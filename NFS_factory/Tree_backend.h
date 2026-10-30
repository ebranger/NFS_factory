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

	//(P mod N_i)^(2^MAX_SQUARINGS_IN_CHECK) mod N_i for every leaf, which Leaf_remainder then gives. The leaves are unchanged.
	virtual void Leaf_remainders() = 0;
	virtual mpz_class* Leaf_remainder(int index) = 0;
};

//A tree backend with the given number of levels (at least 2), for leaves of up to max_bits_in_input bits (to set up the memory).
Tree_backend* create_tree_backend(tree_backend_kind kind, int levels, int max_bits_in_input);

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

//The tree computations on the CPU, with the products of a Tree_products.
class Cpu_tree_backend : public Tree_backend
{
public:
	Cpu_tree_backend(Tree_products* products, int levels, int max_bits_in_input);
	~Cpu_tree_backend();

	mpz_class* Leaf(int index);
	mp_bitcnt_t Build_product_tree();
	void Scaled_root(mpz_srcptr prime_product);
	void Remainder_tree();
	void Leaf_remainders();
	mpz_class* Leaf_remainder(int index);

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
	std::vector<mpz_class> remainder; //of the leaves
	std::vector<mp_bitcnt_t> node_bits; //sizes of the product tree nodes, which set the scale of the fractions.
	mpz_class shifted, inverse, root_remainder, quotient, temp, part, pair_result[2]; //temporaries.
};
