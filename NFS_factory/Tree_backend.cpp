//The tree computations of the batch smoothness check on the CPU, with MPIR's products or with the integer NTT of int_ntt.
//See Bernstein, "How to find smooth parts of integers" and "Scaled remainder trees", for the algorithm.

#include "Tree_backend.h"
#include <stdlib.h>

#define WRAPAROUND_MIN_LIMBS 512 //Use wraparound products in the remainder tree for nodes of at least this many limbs. Tuned on an i7-14700F.
#define TREE_NTT_MIN_LIMBS 512 //The NTT backend multiplies with the NTT when the product has at least this many limbs, and with MPIR below
//(where MPIR is faster, and the exact product is the same).

//An approximation of B^(2n) / a, where a has n limbs and its top bit set, as B^n plus n limbs. It is within a few units of
//floor((B^(2n) - 1) / a), but can be on either side of it.
//One Newton step from the inverse v_h of the top h limbs of a, about half of them: with D = a v_h - B^(n + h), which is less
//than a few B^n in size, the inverse is v_h B^(n - h) - v_h D / B^(2h). Only D modulo 2^L - 1 with L a little more than n limbs
//is needed, and the error of v_h is squared, so the result is as accurate as the truncations of the last step allow.
//D is exact, whatever L the product modulo 2^L - 1 has, so the result does not depend on how the products are computed.
void Tree_products::Newton_inverse(mpz_ptr result, mpz_srcptr a, mp_size_t a_limbs)
{
	if (a_limbs <= INVERSE_EXACT_LIMBS)
	{
		if (result->_mp_alloc < a_limbs + 1)
		{
			mpz_realloc2(result, (a_limbs + 1) * GMP_LIMB_BITS);
		}
		mpn_invert(result->_mp_d, a->_mp_d, a_limbs);
		result->_mp_d[a_limbs] = 1;
		result->_mp_size = (int)(a_limbs + 1);
		return;
	}

	mp_size_t h = a_limbs / 2 + 2;
	mpz_class a_h, v_h, d, power;

	mpz_tdiv_q_2exp(a_h.get_mpz_t(), a, GMP_LIMB_BITS * (a_limbs - h));
	Newton_inverse(v_h.get_mpz_t(), a_h.get_mpz_t(), h);

	//D modulo 2^L - 1, from -2^(L - 1) to 2^(L - 1).
	mp_bitcnt_t mod_bits = Product_mod(d.get_mpz_t(), a, v_h.get_mpz_t(), GMP_LIMB_BITS * (a_limbs + 2));
	mpz_setbit(power.get_mpz_t(), (GMP_LIMB_BITS * (a_limbs + h)) % mod_bits); //B^(n + h) modulo 2^L - 1
	d -= power;
	mpz_set_ui(power.get_mpz_t(), 0);
	mpz_setbit(power.get_mpz_t(), mod_bits);
	power -= 1;
	if (mpz_sgn(d.get_mpz_t()) < 0)
	{
		d += power;
	}
	if (mpz_sizeinbase(d.get_mpz_t(), 2) >= mod_bits)
	{
		d -= power;
	}

	//v_h D / B^(2h), with only the top part of D, which is enough for an error of less than one.
	mpz_fdiv_q_2exp(d.get_mpz_t(), d.get_mpz_t(), GMP_LIMB_BITS * (h - 1));
	Multiply(d.get_mpz_t(), d.get_mpz_t(), v_h.get_mpz_t());
	mpz_fdiv_q_2exp(d.get_mpz_t(), d.get_mpz_t(), GMP_LIMB_BITS * (h + 1));

	mpz_mul_2exp(result, v_h.get_mpz_t(), GMP_LIMB_BITS * (a_limbs - h));
	mpz_sub(result, result, d.get_mpz_t());
}

void Tree_products::Inverse(mpz_ptr result, mpz_srcptr a)
{
	mp_size_t a_limbs = abs(a->_mp_size);
	Newton_inverse(result, a, a_limbs);
	if (a_limbs > INVERSE_EXACT_LIMBS)
	{
		mpz_sub_ui(result, result, INVERSE_MARGIN);
	}
}

//MPIR's products, and its FFT for the products modulo 2^L - 1.
class Mpir_tree_products : public Tree_products
{
public:
	void Multiply(mpz_ptr result, mpz_srcptr a, mpz_srcptr b)
	{
		mpz_mul(result, a, b);
	}

	mp_bitcnt_t Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L)
	{
		return cyclic.Product_mod(result, y, x, min_L);
	}

	void Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision)
	{
		cyclic.Middle_product_pair(result0, result1, y, x0, x1, shift0, shift1, precision);
	}

private:
	Cyclic_product cyclic;
};

//The products of the integer NTT. Small products, whose exact result is the same, are left to MPIR.
class Ntt_tree_products : public Tree_products
{
public:
	void Multiply(mpz_ptr result, mpz_srcptr a, mpz_srcptr b)
	{
		mp_size_t a_limbs = abs(a->_mp_size), b_limbs = abs(b->_mp_size);
		if (a_limbs + b_limbs >= TREE_NTT_MIN_LIMBS && a_limbs >= 2 && b_limbs >= 2)
		{
			ntt.Multiply(result, a, b);
		}
		else
		{
			mpz_mul(result, a, b);
		}
	}

	mp_bitcnt_t Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L)
	{
		return ntt.Product_mod(result, y, x, min_L);
	}

	void Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision)
	{
		ntt.Middle_product_pair(result0, result1, y, x0, x1, shift0, shift1, precision);
	}

private:
	Ntt_product ntt;
};

//The NTT backend: the tree computations of Cpu_tree_backend with the products of the integer NTT, with one plan for each depth
//of the tree (tree_level_plan), and the forward transforms of the nodes from the product tree kept for the middle products of the
//remainder tree. A pair of nodes then costs 3 transforms in the remainder tree instead of 5. The GPU backend computes the same
//numbers in every node.
class Ntt_tree_backend : public Cpu_tree_backend
{
public:
	Ntt_tree_backend(int levels, int max_bits_in_input);

protected:
	void Product_level(int depth);
	void Remainder_level(int depth);

private:
	Ntt_product ntt;
	std::vector<std::vector<uint32_t>> transform; //the forward transforms of the nodes, where kept
	std::vector<Ntt_plan> plan; //of each depth, for the products of its pairs
	std::vector<char> uses_ntt; //of each depth
	std::vector<uint32_t> y_transform, work, work2;
};

#ifdef NFS_GPU
Tree_backend* create_gpu_tree_backend(int levels, int max_bits_in_input); //Gpu_tree_backend.cu
#endif

Tree_backend* create_tree_backend(tree_backend_kind kind, int levels, int max_bits_in_input)
{
#ifdef NFS_GPU
	if (kind == TREE_BACKEND_GPU)
	{
		return create_gpu_tree_backend(levels, max_bits_in_input);
	}
#endif
	if (kind == TREE_BACKEND_NTT)
	{
		return new Ntt_tree_backend(levels, max_bits_in_input);
	}
	return new Cpu_tree_backend(new Mpir_tree_products, levels, max_bits_in_input);
}

//The tree takes ownership of products.
Cpu_tree_backend::Cpu_tree_backend(Tree_products* products, int levels, int max_bits_in_input)
{
	this->products = products;
	this->levels = levels;
	num_leaves = 1 << (levels - 1);
	first_leaf = num_leaves - 1;

	//A guard bit for each level of truncations, and two more.
	guard = levels + 2;

	//Allocate enough memory for every node, so that the numbers do not have to grow during the run.
	node.resize(2 * num_leaves - 1);
	node_bits.resize(2 * num_leaves - 1);
	cofactor.resize(num_leaves);
	mp_bitcnt_t bits = max_bits_in_input;
	for (int depth = levels - 1; depth >= 0; depth--)
	{
		for (int i = (1 << depth) - 1; i <= 2 * ((1 << depth) - 1); i++)
		{
			mpz_realloc2(node[i].get_mpz_t(), bits);
			node[i] = 1;
		}
		bits *= 2;
	}
}

Cpu_tree_backend::~Cpu_tree_backend()
{
	delete products;
}

mpz_class* Cpu_tree_backend::Leaf(int index)
{
	return &node[first_leaf + index];
}

mpz_class* Cpu_tree_backend::Leaf_cofactor(int index)
{
	return &cofactor[index];
}

mp_bitcnt_t Cpu_tree_backend::Build_product_tree()
{
	for (int i = first_leaf; i < first_leaf + num_leaves; i++)
	{
		node_bits[i] = mpz_sizeinbase(node[i].get_mpz_t(), 2);
	}
	for (int depth = levels - 2; depth >= 0; depth--)
	{
		Product_level(depth);
	}
	return node_bits[0];
}

void Cpu_tree_backend::Product_level(int depth)
{
	for (int i = (1 << depth) - 1; i <= 2 * ((1 << depth) - 1); i++)
	{
		Product_node(i);
		node_bits[i] = mpz_sizeinbase(node[i].get_mpz_t(), 2);
	}
}

void Cpu_tree_backend::Product_node(int i)
{
	products->Multiply(node[i].get_mpz_t(), node[2 * i + 1].get_mpz_t(), node[2 * i + 2].get_mpz_t());
}

//result = a modulo 2^L - 1, from 0 to 2^L - 2. temp is used as scratch space.
static void mod_2expm1(mpz_ptr result, mpz_srcptr a, mp_bitcnt_t L, mpz_ptr temp)
{
	mpz_set(result, a);
	while (mpz_sizeinbase(result, 2) > L)
	{
		mpz_tdiv_q_2exp(temp, result, L);
		mpz_tdiv_r_2exp(result, result, L);
		mpz_add(result, result, temp);
	}
	mpz_add_ui(result, result, 1);
	if (mpz_sizeinbase(result, 2) > L)
	{
		mpz_set_ui(result, 0);
	}
	else
	{
		mpz_sub_ui(result, result, 1);
	}
}

//result = the bits of a from position up to position + count, as a number.
static void bits_of(mpz_ptr result, mpz_srcptr a, mp_bitcnt_t position, mp_bitcnt_t count)
{
	mp_size_t a_limbs = abs(a->_mp_size);
	mp_size_t first = (mp_size_t)(position / GMP_LIMB_BITS);
	if (count == 0 || first >= a_limbs)
	{
		mpz_set_ui(result, 0);
		return;
	}

	mp_size_t out_limbs = (mp_size_t)((count + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS);
	mp_size_t source_limbs = (out_limbs + 1 < a_limbs - first) ? out_limbs + 1 : a_limbs - first;
	if (result->_mp_alloc < source_limbs)
	{
		mpz_realloc2(result, source_limbs * GMP_LIMB_BITS);
	}

	if (position % GMP_LIMB_BITS != 0)
	{
		mpn_rshift(result->_mp_d, a->_mp_d + first, source_limbs, (unsigned int)(position % GMP_LIMB_BITS));
	}
	else
	{
		mpn_copyi(result->_mp_d, a->_mp_d + first, source_limbs);
	}

	out_limbs = (out_limbs < source_limbs) ? out_limbs : source_limbs;
	if (count % GMP_LIMB_BITS != 0 && out_limbs == (mp_size_t)((count + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS))
	{
		result->_mp_d[out_limbs - 1] &= (((mp_limb_t)1) << (count % GMP_LIMB_BITS)) - 1;
	}
	while (out_limbs > 0 && result->_mp_d[out_limbs - 1] == 0)
	{
		out_limbs--;
	}
	result->_mp_size = (int)out_limbs;
}

//Replace the root of the product tree X by the root of the scaled remainder tree, floor((P mod X) * 2^precision / X), with
//precision = bits of X + guard bits, where P is the prime product. One inverse of X is computed and used both for P mod X
//(Barrett reduction) and for the scaling, which is faster than two divisions that each compute an inverse. The result can be one
//too small, which the guard bits absorb, just as they absorb the truncation errors further down in the tree.
void Cpu_tree_backend::Scaled_root(mpz_srcptr prime_product)
{
	mpz_ptr root = node[0].get_mpz_t();
	mp_bitcnt_t precision = node_bits[0] + guard;
	mp_bitcnt_t bits = mpz_sizeinbase(root, 2);
	mp_size_t limbs = abs(root->_mp_size);
	mp_bitcnt_t normalize = limbs * GMP_LIMB_BITS - bits;

	//V is the inverse of the root shifted to fill limbs + 1 limbs, root * V < 2^e, and at most INVERSE_MAX_ERROR less than
	//2^e / root, with e = 128 * limbs + 64 - normalize. The extra limb makes the scaling below accurate to less than one unit.
	mpz_mul_2exp(shifted.get_mpz_t(), root, normalize + GMP_LIMB_BITS);
	mpz_ptr v = inverse.get_mpz_t();
	products->Inverse(v, shifted.get_mpz_t());
	mp_bitcnt_t e = 128 * limbs + GMP_LIMB_BITS - normalize;

	//r = P mod X. Barrett reduction works for numbers less than 2^(2 * bits), which the choice of tree size gives. A larger prime
	//product, as in the second stage when its tree is not full, is reduced one part at a time from the top: the remainder so far,
	//followed by the next bits - 1 bits of P, is less than 2^(2 * bits - 1).
	mpz_ptr r = root_remainder.get_mpz_t();
	mp_bitcnt_t product_size = mpz_sizeinbase(prime_product, 2);

	if (product_size <= 2 * bits)
	{
		Barrett_reduce(r, prime_product, root, e, bits);
	}
	else
	{
		mp_bitcnt_t position = product_size;
		mpz_set_ui(r, 0);
		while (position > 0)
		{
			mp_bitcnt_t size = (position > bits - 1) ? bits - 1 : position;
			position -= size;
			bits_of(part.get_mpz_t(), prime_product, position, size);
			mpz_mul_2exp(r, r, size);
			mpz_add(r, r, part.get_mpz_t());
			Barrett_reduce(r, r, root, e, bits);
		}
	}

	//floor(r * 2^precision / X) = floor(r * V / 2^(e - precision)), or one less.
	products->Multiply(root, r, v);
	mpz_tdiv_q_2exp(root, root, e - precision);
}

//r = a mod x, for a less than 2^(2 * bits), with x of bits bits and its inverse V from Scaled_root (inverse), so that x * V < 2^e.
//r can be a.
void Cpu_tree_backend::Barrett_reduce(mpz_ptr r, mpz_srcptr a, mpz_srcptr x, mp_bitcnt_t e, mp_bitcnt_t bits)
{
	mp_bitcnt_t a_size = mpz_sizeinbase(a, 2);
	if (a_size < bits)
	{
		mpz_set(r, a);
		return;
	}

	//mu = floor(2^(2 * bits) / x), or up to two less, without the bits that do not matter for the quotient. The quotient has
	//at most quotient_bits bits, and the bits of mu more than 64 below that change the estimate by less than one.
	//The quotient estimate is then at most 5 too small.
	mpz_ptr q = quotient.get_mpz_t();
	mpz_ptr t = temp.get_mpz_t();
	mp_bitcnt_t quotient_bits = a_size - bits + 1;
	mp_bitcnt_t cut = (bits > quotient_bits + GMP_LIMB_BITS) ? bits - quotient_bits - GMP_LIMB_BITS : 0;
	mpz_tdiv_q_2exp(q, inverse.get_mpz_t(), e - 2 * bits + cut);
	mpz_tdiv_q_2exp(t, a, bits - 1);
	products->Multiply(q, q, t);
	mpz_tdiv_q_2exp(q, q, bits + 1 - cut);

	//r = a - q x is less than 6 x, so it can be computed modulo 2^L - 1 > 6 x, which only needs q x modulo 2^L - 1.
	mp_bitcnt_t L = products->Product_mod(q, q, x, bits + 3);
	mod_2expm1(t, a, L, part.get_mpz_t());
	mpz_sub(r, t, q);
	if (mpz_sgn(r) < 0)
	{
		mpz_set_ui(q, 0);
		mpz_setbit(q, L);
		mpz_sub_ui(q, q, 1);
		mpz_add(r, r, q);
	}
	while (mpz_cmp(r, x) >= 0)
	{
		mpz_sub(r, r, x);
	}
}

//result = ((y * x) mod 2^precision) >> shift, exactly.
void Cpu_tree_backend::Truncated_product(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t shift, mp_bitcnt_t precision)
{
	mpz_mul(result, y, x);
	mpz_tdiv_r_2exp(result, result, precision);
	mpz_tdiv_q_2exp(result, result, shift);
}

//A node of N bits whose parent has a fraction y of parent_bits + guard bits (the fraction times 2^(parent_bits + guard)) gets
//the fraction y * sibling, which is scaled by 2^(parent_bits + guard), and keeps its top N + guard bits: the bits from
//parent_bits - N up to parent_bits + guard. The bits above are a whole number, which does not matter for a fraction modulo 1.
void Cpu_tree_backend::Remainder_tree()
{
	for (int depth = 1; depth <= levels - 2; depth++)
	{
		Remainder_level(depth);
	}
}

void Cpu_tree_backend::Remainder_level(int depth)
{
	for (int i = (1 << depth) - 1; i <= 2 * ((1 << depth) - 1); i += 2)
	{
		Remainder_pair(i);
	}
}

void Cpu_tree_backend::Remainder_pair(int i)
{
	mpz_srcptr parent = node[(i - 1) / 2].get_mpz_t();
	mp_bitcnt_t parent_bits = node_bits[(i - 1) / 2];
	mp_bitcnt_t precision = parent_bits + guard;
	if (precision / GMP_LIMB_BITS >= WRAPAROUND_MIN_LIMBS)
	{
		//Large numbers: only the middle part of the products is needed, which is faster to get with a wraparound product.
		products->Middle_product_pair(pair_result[0].get_mpz_t(), pair_result[1].get_mpz_t(), parent, node[i + 1].get_mpz_t(), node[i].get_mpz_t(),
			parent_bits - node_bits[i], parent_bits - node_bits[i + 1], precision);
	}
	else
	{
		Truncated_product(pair_result[0].get_mpz_t(), parent, node[i + 1].get_mpz_t(), parent_bits - node_bits[i], precision);
		Truncated_product(pair_result[1].get_mpz_t(), parent, node[i].get_mpz_t(), parent_bits - node_bits[i + 1], precision);
	}
	mpz_swap(node[i].get_mpz_t(), pair_result[0].get_mpz_t());
	mpz_swap(node[i + 1].get_mpz_t(), pair_result[1].get_mpz_t());
}

//At a leaf N_i, the fraction y (scaled by 2^(N_i bits + guard)) times N_i is P mod N_i, scaled by 2^guard. The truncations only
//made y smaller, by less than the guard bits allow, so rounding up gives P mod N_i. It is then squared modulo N_i, and the
//cofactor is N_i / gcd(N_i, Q).
void Cpu_tree_backend::Leaf_cofactors()
{
#pragma omp parallel for schedule(dynamic, 256)
	for (int i = first_leaf; i < first_leaf + num_leaves; i++)
	{
		const int parent = (i - 1) / 2;
		const int sibling = (i % 2 == 1) ? i + 1 : i - 1;
		mpz_ptr q = cofactor[i - first_leaf].get_mpz_t();
		mpz_srcptr n = node[i].get_mpz_t();

		Truncated_product(q, node[parent].get_mpz_t(), node[sibling].get_mpz_t(), node_bits[parent] - node_bits[i], node_bits[parent] + guard);
		mpz_mul(q, q, n);
		mpz_tdiv_q_2exp(q, q, node_bits[i]);
		mpz_add_ui(q, q, (1UL << guard) - 1);
		mpz_tdiv_q_2exp(q, q, guard);

		for (int j = 0; j < MAX_SQUARINGS_IN_CHECK; j++)
		{
			mpz_mul(q, q, q);
			mpz_tdiv_r(q, q, n);
		}

		mpz_class common;
		mpz_gcd(common.get_mpz_t(), n, q);
		mpz_divexact(q, n, common.get_mpz_t());
	}
}

bool tree_level_plan(Ntt_plan& plan, const mp_bitcnt_t* bits, size_t count, mp_bitcnt_t guard)
{
	mp_bitcnt_t max_pair = 0, max_node = 0;
	for (size_t i = 0; i + 1 < count; i += 2)
	{
		max_pair = (bits[i] + bits[i + 1] > max_pair) ? bits[i] + bits[i + 1] : max_pair;
		max_node = (bits[i] > max_node) ? bits[i] : max_node;
		max_node = (bits[i + 1] > max_node) ? bits[i + 1] : max_node;
	}
	if (max_node < TREE_NTT_MIN_CHILD_BITS)
	{
		return false;
	}
	return ntt_plan_cyclic(plan, max_pair + guard + GMP_LIMB_BITS, max_pair + guard, max_node);
}

Ntt_tree_backend::Ntt_tree_backend(int levels, int max_bits_in_input) : Cpu_tree_backend(new Ntt_tree_products, levels, max_bits_in_input)
{
	transform.resize(node.size());
	plan.resize(levels);
	uses_ntt.assign(levels, 0);
}

void Ntt_tree_backend::Product_level(int depth)
{
	const int child_depth = depth + 1;
	const int first_child = (1 << child_depth) - 1, count = 1 << child_depth;
	uses_ntt[child_depth] = tree_level_plan(plan[child_depth], &node_bits[first_child], count, guard);
	if (!uses_ntt[child_depth])
	{
		Cpu_tree_backend::Product_level(depth);
		return;
	}
	const Ntt_plan& p = plan[child_depth];
	const bool keep = child_depth <= levels - 2; //The leaves are not needed in the remainder tree.
	for (int i = first_child; i < first_child + count; i += 2)
	{
		std::vector<uint32_t>& left = keep ? transform[i] : work2;
		std::vector<uint32_t>& right = keep ? transform[i + 1] : y_transform;
		ntt.Forward(left, p, node[i].get_mpz_t());
		ntt.Forward(right, p, node[i + 1].get_mpz_t());
		work = left;
		ntt.Pointwise(work, right, p);
		const int parent = (i - 1) / 2;
		ntt.Product_from(node[parent].get_mpz_t(), work, p);
		node_bits[parent] = mpz_sizeinbase(node[parent].get_mpz_t(), 2);
	}
}

void Ntt_tree_backend::Remainder_level(int depth)
{
	const int first = (1 << depth) - 1, count = 1 << depth;
	for (int i = first; i < first + count; i += 2)
	{
		const int parent = (i - 1) / 2;
		const mp_bitcnt_t precision = node_bits[parent] + guard;
		if (!uses_ntt[depth])
		{
			Truncated_product(pair_result[0].get_mpz_t(), node[parent].get_mpz_t(), node[i + 1].get_mpz_t(), node_bits[parent] - node_bits[i], precision);
			Truncated_product(pair_result[1].get_mpz_t(), node[parent].get_mpz_t(), node[i].get_mpz_t(), node_bits[parent] - node_bits[i + 1], precision);
		}
		else if (mpz_sgn(node[parent].get_mpz_t()) == 0)
		{
			pair_result[0] = 0;
			pair_result[1] = 0;
		}
		else
		{
			const Ntt_plan& p = plan[depth];
			ntt.Forward(y_transform, p, node[parent].get_mpz_t());
			work = transform[i + 1];
			ntt.Pointwise(work, y_transform, p);
			ntt.Middle_from(pair_result[0].get_mpz_t(), work, p, node_bits[parent] - node_bits[i], precision);
			work = transform[i];
			ntt.Pointwise(work, y_transform, p);
			ntt.Middle_from(pair_result[1].get_mpz_t(), work, p, node_bits[parent] - node_bits[i + 1], precision);
		}
		mpz_swap(node[i].get_mpz_t(), pair_result[0].get_mpz_t());
		mpz_swap(node[i + 1].get_mpz_t(), pair_result[1].get_mpz_t());
	}
}
