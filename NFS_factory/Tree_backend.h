#pragma once
#include "mpirxx.h"
#include <vector>
#include "Fft_product.h"
#include "int_ntt.h"

#define MAX_SQUARINGS_IN_CHECK 1 //The leaf remainders are squared this many times modulo their leaf, so that prime powers up to 2^this are found.
//Currently: do only 1 squaring, it misses a few very rare relations but increases speed slightly.

//The implementations of the tree computations.
enum tree_backend_kind
{
	TREE_BACKEND_MPIR, //MPIR's products, with the products modulo 2^L - 1 of Fft_product.
	TREE_BACKEND_NTT,  //The integer NTT of int_ntt on the CPU, the reference for the GPU version: the same numbers in every node.
	TREE_BACKEND_GPU   //The integer NTT on the GPU (only in the GPU build).
};

#define TREE_NTT_MIN_CHILD_BITS 1024 //The NTT and GPU backends multiply the pairs of a depth of the tree with the NTT if the largest
//node there has at least this many bits, and exactly (schoolbook on the GPU) below.

//The plan of the NTT and GPU backends for the count nodes (pairs) of a depth of the tree with the given sizes. It is valid for the
//product of each pair, and for the middle products of the remainder tree there, which need L at least
//max(precision + bits of x - shift, precision) + 64 (see Ntt_product::Middle_from), with precision = parent bits + guard and
//shift = parent bits - node bits: at most the bits of the pair + guard + 64, as the parent has at most the bits of the pair.
//False if the depth uses exact products.
bool tree_level_plan(Ntt_plan& plan, const mp_bitcnt_t* bits, size_t count, mp_bitcnt_t guard);

//The checks of the cofactors of the leaves that a backend can do itself, as Batch_smooth::Classify_cofactor, with its limits:
//cofactors of more than max_cofactor_bits bits are 3LP candidates if composite, of at most three_lp_max_bits bits and at least
//three_lp_min (two limbs); those of more than two_lp_min_bits bits are split into two large primes of at most one_lp_bits bits;
//and those of at most one_lp_bits bits are one large prime if prime and more than max_batch_prime.
struct Leaf_check
{
	int max_cofactor_bits, three_lp_max_bits, two_lp_min_bits, one_lp_bits;
	unsigned long long three_lp_min[2];
	unsigned long long max_batch_prime;
};

//The result of the check of a leaf: the first four as Batch_smooth::cofactor_result, and the cofactors that the CPU checks.
enum leaf_class { LEAF_NOT_SMOOTH, LEAF_SMOOTH, LEAF_CANDIDATE, LEAF_TINYQS, LEAF_CHECK_ON_CPU };

//The batch smoothness check of the numbers N_i of a batch against a prime product P, with Bernstein's scaled remainder tree.
//The product tree has the N_i as leaves and their product X at the root. Each node of the remainder tree then holds the
//fraction (P mod N) / N of its node value N, as a fixed-point number: the root from P mod X, and each other node from its parent's
//fraction times its sibling, since P / N = (P / parent) * sibling. Only the bits that matter for the leaves are kept, and the
//truncation errors are absorbed by guard bits. At a leaf, the fraction times N_i, rounded, is P mod N_i.
//The steps work on whole tree levels, so that a GPU backend can do each level with batched launches.
class Tree_backend
{
public:
	virtual ~Tree_backend() {}

	//Leaf i of the tree, from 0 to 2^(levels - 1) - 1. The caller sets the leaves before Build_product_tree. Unused leaves must be 1.
	virtual mpz_class* Leaf(int index) = 0;

	//The product tree of the leaves. Returns the size of the root in bits.
	virtual mp_bitcnt_t Build_product_tree() = 0;

	//The root of the remainder tree from the prime product, and the remainder tree down to the parents of the leaves.
	virtual void Scaled_root(mpz_srcptr prime_product) = 0;
	virtual void Remainder_tree() = 0;

	//The cofactor N_i / gcd(N_i, (P mod N_i)^(2^MAX_SQUARINGS_IN_CHECK)) of every leaf, which Leaf_cofactor then gives: the part of
	//N_i without the primes of P (that occur to at most the power 2^MAX_SQUARINGS_IN_CHECK), 1 if N_i is smooth. The leaves are
	//unchanged.
	virtual void Leaf_cofactors() = 0;
	virtual mpz_class* Leaf_cofactor(int index) = 0;

	//Have Leaf_cofactors check the cofactors too, if the backend can; false if it cannot. After Leaf_cofactors, Leaf_classes then
	//gives the class of each leaf (or NULL if this batch was not checked), and Leaf_large_primes its two large primes for
	//LEAF_SMOOTH (0 for none). Leaf_cofactor is then only set for the leaves that are not LEAF_SMOOTH or LEAF_NOT_SMOOTH.
	virtual bool Set_leaf_check(const Leaf_check& check) { return false; }
	virtual const unsigned char* Leaf_classes() { return NULL; }
	virtual const unsigned long long* Leaf_large_primes() { return NULL; }

	//Have Leaf_cofactors also find the prime factors of the smooth leaves for the output, if the backend can (the GPU, in the
	//background); false if it cannot. The smooth leaves are those of LEAF_SMOOTH if the backend checks the cofactors, and
	//otherwise those with cofactor 1. Leaf_factors_ready waits for the factors, and returns false if there are none for this
	//batch. Leaf_factors then gives the factors of leaf i without its large primes, in increasing order, and returns their number,
	//or -1 if the leaf was not factored. Leaf_factors is thread safe.
	virtual bool Set_leaf_factoring() { return false; }
	virtual bool Leaf_factors_ready() { return false; }
	virtual int Leaf_factors(int index, const unsigned int** factors) { return -1; }
};

//A tree backend with the given number of levels (at least 2), for leaves of up to max_bits_in_input bits (to set up the memory).
Tree_backend* create_tree_backend(tree_backend_kind kind, int levels, int max_bits_in_input);

//The default number of batches in progress at once for a backend (see Batch_smooth), each with a tree of this size: 1 on the CPU,
//and on the GPU as many as fit in its free memory, up to 6 (1 if not even two fit). Also no more than fit in half of the free
//memory of the computer.
int default_pipeline(tree_backend_kind kind, int levels, int max_bits_in_input);

#ifdef NFS_GPU
//Use GPU number device (from 0) for the gpu backend. Returns false, after saying why, if there is no such CUDA GPU, or if it is
//older than compute capability 7.5. Call it before anything else uses the GPU; otherwise GPU 0 is used.
bool gpu_init(int device);
#endif

//The large products of the tree computations on the CPU.
class Tree_products
{
public:
	virtual ~Tree_products() {}

	//result = a * b.
	virtual void Multiply(mpz_ptr result, mpz_srcptr a, mpz_srcptr b) = 0;

	//result = y * x modulo 2^L - 1, from 0 to 2^L - 2, with L at least min_L. Both numbers must be less than 2^L. Returns L.
	virtual mp_bitcnt_t Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L) = 0;

	//result0 = ((y * x0) mod 2^precision) >> shift0, and result1 the same from x1 and shift1, or one less (modulo 2^(precision - shift)).
	//y must be less than 2^precision, and x0 and x1 less than 2^(shift + 1).
	virtual void Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision) = 0;

	//An approximation of B^(2n) / a, where a has n limbs and its top bit set, with B = 2^64. The result is at most
	//INVERSE_MAX_ERROR less than floor((B^(2n) - 1) / a), and never more. It does not depend on how the products are computed.
	void Inverse(mpz_ptr result, mpz_srcptr a);

private:
	void Newton_inverse(mpz_ptr result, mpz_srcptr a, mp_size_t a_limbs);
};

#define INVERSE_MAX_ERROR (1 << 16)
#define INVERSE_EXACT_LIMBS 2000 //Inverses of up to this many limbs are computed exactly by MPIR.
#define INVERSE_MARGIN (1 << 12) //Subtracted from the Newton inverse, which can be a few units too large.

//The tree computations on the CPU, with the products of a Tree_products.
class Cpu_tree_backend : public Tree_backend
{
public:
	//With preallocate_nodes, the memory of all nodes is allocated up front; otherwise only that of the leaves (the GPU backend
	//keeps the other nodes on the GPU).
	Cpu_tree_backend(Tree_products* products, int levels, int max_bits_in_input, bool preallocate_nodes = true);
	~Cpu_tree_backend();

	mpz_class* Leaf(int index);
	mp_bitcnt_t Build_product_tree();
	void Scaled_root(mpz_srcptr prime_product);
	void Remainder_tree();
	void Leaf_cofactors();
	mpz_class* Leaf_cofactor(int index);

protected:
	//The products of the nodes at a depth from their children, and the fractions of the nodes at a depth from their parents'.
	//Product_level sets the sizes in node_bits.
	virtual void Product_level(int depth);
	virtual void Remainder_level(int depth);
	void Product_node(int i);
	void Remainder_pair(int i);
	void Barrett_reduce(mpz_ptr r, mpz_srcptr a, mpz_srcptr x, mp_bitcnt_t e, mp_bitcnt_t bits);
	void Truncated_product(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t shift, mp_bitcnt_t precision);

	Tree_products* products;
	int levels;
	int first_leaf; //the nodes are numbered from the root (0), with the children of node i at 2i + 1 and 2i + 2.
	int num_leaves;
	mp_bitcnt_t guard; //guard bits of the fixed-point fractions.
	std::vector<mpz_class> node;
	std::vector<mpz_class> cofactor; //of the leaves
	std::vector<mp_bitcnt_t> node_bits; //sizes of the product tree nodes, which set the scale of the fractions.
	mpz_class shifted, inverse, root_remainder, quotient, temp, part, pair_result[2]; //temporaries.
};
