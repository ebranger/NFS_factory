//The tree computations on the GPU (TREE_BACKEND_GPU, in the GPU build only). Each depth of the product and remainder trees is
//done with batched kernels for all its nodes: the NTT products of int_ntt with one plan per depth (tree_level_plan), or schoolbook
//products for small nodes, exactly as Ntt_tree_backend does on the CPU, so that every node gets the same number. The forward
//transforms of the nodes from the product tree are kept on the GPU for the middle products of the remainder tree.
//The root of the remainder tree is computed on the GPU too (Gpu_root), with the same numbers as Cpu_tree_backend::Scaled_root,
//and the cofactors of the leaves, one thread per leaf (Gpu_ntt::Leaf_cofactors), of which only the cofactors are copied back.
//Leaves of more than 8 limbs are done on the CPU, by the code of Cpu_tree_backend, from the fractions of their parents.

#include "Tree_backend.h"
#include "ntt_gpu.h"
#include "Parallel.h"
#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <algorithm>

//Products on the GPU one at a time, with copies to and from the GPU, as the Tree_products of the code of Cpu_tree_backend that
//the GPU tree uses (the GPU tree computes its root with Gpu_root instead). Small products are left to MPIR.
class Gpu_tree_products : public Tree_products
{
public:
	Gpu_tree_products()
	{
		in_limbs = 0;
		out_limbs = 0;
		transform_words = 0;
		limbs[0] = limbs[1] = NULL;
		sizes = NULL;
		X = Y = NULL;
		CUDA_CHECK(cudaMalloc(&sizes, 4 * sizeof(uint32_t)));
	}

	~Gpu_tree_products()
	{
		cudaFree(limbs[0]);
		cudaFree(limbs[1]);
		cudaFree(sizes);
		cudaFree(X);
		cudaFree(Y);
	}

	void Multiply(mpz_ptr result, mpz_srcptr a, mpz_srcptr b)
	{
		mp_size_t a_limbs = abs(a->_mp_size), b_limbs = abs(b->_mp_size);
		Ntt_plan plan;
		if (a_limbs + b_limbs < 512 || a_limbs < 2 || b_limbs < 2 || !ntt_plan_product(plan, mpz_sizeinbase(a, 2), mpz_sizeinbase(b, 2)))
		{
			mpz_mul(result, a, b);
			return;
		}
		const bool negative = (a->_mp_size < 0) != (b->_mp_size < 0);
		Product(result, a, b, plan, false);
		if (negative)
		{
			mpz_neg(result, result);
		}
	}

	mp_bitcnt_t Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L)
	{
		Ntt_plan plan;
		if (!ntt_plan_cyclic(plan, min_L, mpz_sizeinbase(y, 2), mpz_sizeinbase(x, 2)))
		{
			return cpu.Product_mod(result, y, x, min_L);
		}
		if (mpz_sgn(y) == 0 || mpz_sgn(x) == 0)
		{
			mpz_set_ui(result, 0);
			return plan.L();
		}
		Product(result, y, x, plan, true);
		return plan.L();
	}

	void Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision)
	{
		cpu.Middle_product_pair(result0, result1, y, x0, x1, shift0, shift1, precision); //Not used by the GPU tree.
	}

	Gpu_ntt gpu;

private:
	//result = a * b, exactly or modulo 2^L - 1, with the plan.
	void Product(mpz_ptr result, mpz_srcptr a, mpz_srcptr b, const Ntt_plan& plan, bool cyclic)
	{
		size_t a_limbs = mpz_size(a), b_limbs = mpz_size(b);
		size_t stride = (a_limbs > b_limbs) ? a_limbs : b_limbs;
		size_t out_stride = Gpu_ntt::Combine_limbs(plan);
		Reserve(2 * stride, out_stride, (size_t)plan.primes * 2 * plan.Length());
		CUDA_CHECK(cudaMemcpy(limbs[0], a->_mp_d, a_limbs * sizeof(uint64_t), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(limbs[0] + stride, b->_mp_d, b_limbs * sizeof(uint64_t), cudaMemcpyHostToDevice));
		uint32_t size[2] = { (uint32_t)a_limbs, (uint32_t)b_limbs };
		CUDA_CHECK(cudaMemcpy(sizes, size, sizeof(size), cudaMemcpyHostToDevice));
		Gpu_level in = { limbs[0], stride, sizes, 2 };
		Gpu_level out = { limbs[1], out_stride, sizes + 2, 1 };
		gpu.Forward(plan, in, X, 0);
		gpu.Pointwise_adjacent(plan, X, Y, 1, 0);
		gpu.Combine(plan, Y, 1, cyclic, out, 0);
		uint32_t n;
		CUDA_CHECK(cudaMemcpy(&n, sizes + 2, sizeof(uint32_t), cudaMemcpyDeviceToHost));
		if (result->_mp_alloc < (int)n)
		{
			mpz_realloc2(result, n * GMP_LIMB_BITS);
		}
		CUDA_CHECK(cudaMemcpy(result->_mp_d, limbs[1], n * sizeof(uint64_t), cudaMemcpyDeviceToHost));
		result->_mp_size = (int)n;
	}

	void Reserve(size_t in, size_t out, size_t words)
	{
		if (in > in_limbs)
		{
			cudaFree(limbs[0]);
			CUDA_CHECK(cudaMalloc(&limbs[0], in * sizeof(uint64_t)));
			in_limbs = in;
		}
		if (out > out_limbs)
		{
			cudaFree(limbs[1]);
			CUDA_CHECK(cudaMalloc(&limbs[1], out * sizeof(uint64_t)));
			out_limbs = out;
		}
		if (words > transform_words)
		{
			cudaFree(X);
			cudaFree(Y);
			CUDA_CHECK(cudaMalloc(&X, words * sizeof(uint32_t)));
			CUDA_CHECK(cudaMalloc(&Y, words / 2 * sizeof(uint32_t)));
			transform_words = words;
		}
	}

	Ntt_product cpu; //for products beyond the longest transform
	uint64_t* limbs[2];
	uint32_t* sizes;
	uint32_t *X, *Y;
	size_t in_limbs, out_limbs, transform_words;
};

//Device memory that grows as needed and is kept.
template <class T> static void reserve(T*& p, size_t& capacity, size_t count)
{
	if (count > capacity)
	{
		cudaFree(p);
		CUDA_CHECK(cudaMalloc(&p, count * sizeof(T)));
		capacity = count;
	}
}

//A number on the GPU: n limbs, an upper bound of its size, with zeros above the value. The memory grows as needed, and its
//contents are lost when it grows.
struct Gpu_number
{
	uint64_t* d;
	size_t n, capacity;

	Gpu_number() : d(NULL), n(0), capacity(0) {}
	~Gpu_number() { cudaFree(d); }
	Gpu_number(const Gpu_number&) = delete;
	Gpu_number& operator=(const Gpu_number&) = delete;

	void Resize(size_t limbs)
	{
		reserve(d, capacity, (limbs > 0) ? limbs : 1);
		n = limbs;
	}

	void Swap(Gpu_number& other)
	{
		std::swap(d, other.d);
		std::swap(n, other.n);
		std::swap(capacity, other.capacity);
	}
};

#define ROOT_MAX_DEPTH 40 //of the Newton iteration, which halves the length down to INVERSE_EXACT_LIMBS

//The root of the scaled remainder tree on the GPU: the steps of Cpu_tree_backend::Scaled_root, Barrett_reduce and
//Tree_products::Inverse on numbers in device memory, with the same numbers in every step, so that the result is the same.
//Only the sizes of numbers are copied to the CPU, to choose the plans, and the top INVERSE_EXACT_LIMBS limbs of the root for
//the start of the Newton iteration.
class Gpu_root
{
public:
	Gpu_root(Gpu_ntt* gpu);
	~Gpu_root();

	//fraction = floor((P mod X) 2^precision / X), or one less, from X on the GPU (x_limbs limbs without leading zeros, bits
	//bits) and the prime product P.
	void Compute(Gpu_number& fraction, const uint64_t* X, size_t x_limbs, mp_bitcnt_t bits, mp_bitcnt_t precision, mpz_srcptr prime_product);

private:
	void Newton_inverse(int depth, const uint64_t* a, size_t n);
	const uint64_t* Floor_shift(Gpu_number& x, size_t limbs, bool negative, size_t& out_limbs);
	void Barrett_reduce(Gpu_number& r, const uint64_t* a, size_t a_limbs, mp_bitcnt_t a_bits, const uint64_t* x, size_t x_limbs, mp_bitcnt_t e, mp_bitcnt_t bits);
	void Mod_2expm1(Gpu_number& result, const uint64_t* a, size_t a_limbs, mp_bitcnt_t L);
	void Multiply(Gpu_number& out, const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs);
	mp_bitcnt_t Product_mod(Gpu_number& out, const uint64_t* y, size_t y_limbs, const uint64_t* x, size_t x_limbs, mp_bitcnt_t min_L);
	bool Transform_pair(const Ntt_plan& plan, const uint64_t* a, const uint64_t* b);
	void Sizes(const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs);
	mp_bitcnt_t Bits(const uint64_t* a, size_t limbs);
	uint32_t Carry();
	void Zero(Gpu_number& x, size_t limbs);
	void To_mpz(mpz_ptr result, const uint64_t* a, size_t limbs);
	void From_mpz(Gpu_number& x, mpz_srcptr a, size_t min_limbs);

	Gpu_ntt* gpu;
	uint32_t* words; //on the GPU: sizes of the operands of a product (limbs, bits) at 0 to 3, of the result at 4, and a carry at 6.
	uint32_t host_words[8];
	uint32_t *X, *Y;
	size_t transform_words;
	Gpu_number level[ROOT_MAX_DEPTH]; //the inverses of the Newton iteration, of the whole root at 0
	Gpu_number shifted, d, magnitude, product, remainder, part, quotient, temp, temp2, prime_product_gpu;
	mpz_srcptr prime_product_host; //of prime_product_gpu
	size_t prime_product_limbs;
	std::vector<mp_limb_t> host_a, host_v;
	mpz_class host_x, host_y, host_z;
};

Gpu_root::Gpu_root(Gpu_ntt* gpu)
{
	this->gpu = gpu;
	CUDA_CHECK(cudaMalloc(&words, 8 * sizeof(uint32_t)));
	X = Y = NULL;
	transform_words = 0;
	prime_product_host = NULL;
	prime_product_limbs = 0;
}

Gpu_root::~Gpu_root()
{
	cudaFree(words);
	cudaFree(X);
	cudaFree(Y);
}

void Gpu_root::Zero(Gpu_number& x, size_t limbs)
{
	x.Resize(limbs);
	CUDA_CHECK(cudaMemsetAsync(x.d, 0, x.n * sizeof(uint64_t), 0));
}

void Gpu_root::To_mpz(mpz_ptr result, const uint64_t* a, size_t limbs)
{
	if (result->_mp_alloc < (int)limbs + 1)
	{
		mpz_realloc2(result, (limbs + 1) * GMP_LIMB_BITS);
	}
	CUDA_CHECK(cudaMemcpy(result->_mp_d, a, limbs * sizeof(uint64_t), cudaMemcpyDeviceToHost));
	while (limbs > 0 && result->_mp_d[limbs - 1] == 0)
	{
		limbs--;
	}
	result->_mp_size = (int)limbs;
}

void Gpu_root::From_mpz(Gpu_number& x, mpz_srcptr a, size_t min_limbs)
{
	const size_t limbs = mpz_size(a);
	Zero(x, (limbs > min_limbs) ? limbs : min_limbs);
	CUDA_CHECK(cudaMemcpy(x.d, a->_mp_d, limbs * sizeof(uint64_t), cudaMemcpyHostToDevice));
}

//The sizes of a and b: limbs and bits at words 0 and 1, and 2 and 3, also in host_words.
void Gpu_root::Sizes(const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs)
{
	gpu->Size(a, a_limbs, words, words + 1, 0);
	gpu->Size(b, b_limbs, words + 2, words + 3, 0);
	CUDA_CHECK(cudaMemcpy(host_words, words, 4 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
}

mp_bitcnt_t Gpu_root::Bits(const uint64_t* a, size_t limbs)
{
	gpu->Size(a, limbs, words, words + 1, 0);
	CUDA_CHECK(cudaMemcpy(host_words, words, 2 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
	return host_words[1];
}

uint32_t Gpu_root::Carry()
{
	CUDA_CHECK(cudaMemcpy(host_words + 6, words + 6, sizeof(uint32_t), cudaMemcpyDeviceToHost));
	return host_words[6];
}

//The pointwise product of the transforms of a and b in X, whose sizes Sizes has set.
bool Gpu_root::Transform_pair(const Ntt_plan& plan, const uint64_t* a, const uint64_t* b)
{
	const size_t words_needed = (size_t)plan.primes * plan.Length();
	if (words_needed > transform_words)
	{
		cudaFree(X);
		cudaFree(Y);
		CUDA_CHECK(cudaMalloc(&X, words_needed * sizeof(uint32_t)));
		CUDA_CHECK(cudaMalloc(&Y, words_needed * sizeof(uint32_t)));
		transform_words = words_needed;
	}
	Gpu_level la = { const_cast<uint64_t*>(a), host_words[0], words, 1 };
	Gpu_level lb = { const_cast<uint64_t*>(b), host_words[2], words + 2, 1 };
	gpu->Forward(plan, la, X, 0);
	gpu->Forward(plan, lb, Y, 0);
	gpu->Pointwise(plan, X, Y, 1, 0);
	return true;
}

//out = a * b. out must not be a or b.
void Gpu_root::Multiply(Gpu_number& out, const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs)
{
	Sizes(a, a_limbs, b, b_limbs);
	Ntt_plan plan;
	if (host_words[0] == 0 || host_words[2] == 0)
	{
		Zero(out, a_limbs + b_limbs);
		return;
	}
	if (!ntt_plan_product(plan, host_words[1], host_words[3]))
	{
		//Beyond the longest transform.
		To_mpz(host_x.get_mpz_t(), a, host_words[0]);
		To_mpz(host_y.get_mpz_t(), b, host_words[2]);
		mpz_mul(host_z.get_mpz_t(), host_x.get_mpz_t(), host_y.get_mpz_t());
		From_mpz(out, host_z.get_mpz_t(), a_limbs + b_limbs);
		return;
	}
	Transform_pair(plan, a, b);
	out.Resize(Gpu_ntt::Combine_limbs(plan));
	Gpu_level o = { out.d, out.n, words + 4, 1 };
	gpu->Combine(plan, X, 1, false, o, 0);
}

//out = y * x modulo 2^L - 1, from 0 to 2^L - 2, with the plan of Gpu_tree_products::Product_mod. Returns L. out must not be x or y.
mp_bitcnt_t Gpu_root::Product_mod(Gpu_number& out, const uint64_t* y, size_t y_limbs, const uint64_t* x, size_t x_limbs, mp_bitcnt_t min_L)
{
	Sizes(y, y_limbs, x, x_limbs);
	Ntt_plan plan;
	if (!ntt_plan_cyclic(plan, min_L, host_words[1], host_words[3]))
	{
		To_mpz(host_x.get_mpz_t(), y, host_words[0]);
		To_mpz(host_y.get_mpz_t(), x, host_words[2]);
		Ntt_product cpu;
		mp_bitcnt_t L = cpu.Product_mod(host_z.get_mpz_t(), host_x.get_mpz_t(), host_y.get_mpz_t(), min_L);
		From_mpz(out, host_z.get_mpz_t(), (size_t)(L / GMP_LIMB_BITS + 2));
		return L;
	}
	if (host_words[0] == 0 || host_words[2] == 0)
	{
		Zero(out, Gpu_ntt::Combine_limbs(plan));
		return plan.L();
	}
	Transform_pair(plan, y, x);
	out.Resize(Gpu_ntt::Combine_limbs(plan));
	Gpu_level o = { out.d, out.n, words + 4, 1 };
	gpu->Combine(plan, X, 1, true, o, 0);
	return plan.L();
}

//The inverse of the n limbs at a (top bit set) as Tree_products::Newton_inverse, into level[depth]: an approximation of
//B^(2n) / a as B^n plus n limbs, in n + 2 limbs.
void Gpu_root::Newton_inverse(int depth, const uint64_t* a, size_t n)
{
	if (depth >= ROOT_MAX_DEPTH)
	{
		fprintf(stderr, "Gpu_root: the Newton iteration is too deep\n");
		exit(3);
	}
	Gpu_number& v = level[depth];
	if (n <= INVERSE_EXACT_LIMBS)
	{
		host_a.resize(n);
		host_v.resize(n + 2);
		CUDA_CHECK(cudaMemcpy(host_a.data(), a, n * sizeof(uint64_t), cudaMemcpyDeviceToHost));
		mpn_invert(host_v.data(), host_a.data(), n);
		host_v[n] = 1;
		host_v[n + 1] = 0;
		v.Resize(n + 2);
		CUDA_CHECK(cudaMemcpy(v.d, host_v.data(), (n + 2) * sizeof(uint64_t), cudaMemcpyHostToDevice));
		return;
	}

	const size_t h = n / 2 + 2;
	Newton_inverse(depth + 1, a + (n - h), h);
	Gpu_number& v_h = level[depth + 1];

	//D modulo 2^L - 1, from -2^(L - 1) to 2^(L - 1): a v_h - B^(n + h), with B^(n + h) modulo 2^L - 1 = 2^k.
	const mp_bitcnt_t L = Product_mod(d, a, n, v_h.d, v_h.n, GMP_LIMB_BITS * (n + 2));
	const mp_bitcnt_t k = (GMP_LIMB_BITS * (n + h)) % L;
	gpu->Add(d.d, d.n, d.d, d.n, NULL, 0, (size_t)(k / 64), 1ULL << (k % 64), true, words + 6, 0);
	if (Carry())
	{
		gpu->Add(d.d, d.n, d.d, d.n, NULL, 0, (size_t)(L / 64), 1ULL << (L % 64), false, NULL, 0);
		gpu->Add(d.d, d.n, d.d, d.n, NULL, 0, 0, 1, true, NULL, 0);
	}
	const bool negative = Bits(d.d, d.n) >= L;
	Gpu_number* D = &d;
	if (negative)
	{
		//|D| = 2^L - 1 - d
		magnitude.Resize(d.n);
		gpu->Shift(magnitude.d, magnitude.n, d.d, d.n, 0, L, true, 0);
		D = &magnitude;
	}

	//v_h D / B^(2h), with only the top part of D, as floor(floor(D / B^(h - 1)) v_h / B^(h + 1)).
	size_t d1_limbs, d3_limbs;
	const uint64_t* d1 = Floor_shift(*D, h - 1, negative, d1_limbs);
	Multiply(product, d1, d1_limbs, v_h.d, v_h.n);
	const uint64_t* d3 = Floor_shift(product, h + 1, negative, d3_limbs);

	//v = v_h B^(n - h) - D v_h / B^(2h)
	Zero(v, n + 2);
	const size_t copy = (v_h.n < h + 2) ? v_h.n : h + 2;
	CUDA_CHECK(cudaMemcpyAsync(v.d + (n - h), v_h.d, copy * sizeof(uint64_t), cudaMemcpyDeviceToDevice, 0));
	gpu->Add(v.d, v.n, v.d, v.n, d3, (d3_limbs < v.n) ? d3_limbs : v.n, 0, 0, !negative, NULL, 0);
}

//The magnitude of floor(D / B^limbs) for D = x or D = -x, as limbs at the returned pointer: x / B^limbs for x, and
//ceil(x / B^limbs) = floor((x - 1) / B^limbs) + 1 for -x (x > 0). x is changed. x must have a zero limb at the top.
const uint64_t* Gpu_root::Floor_shift(Gpu_number& x, size_t limbs, bool negative, size_t& out_limbs)
{
	if (x.n <= limbs + 1)
	{
		//Not with the sizes of the Newton iteration, where x has more limbs.
		To_mpz(host_x.get_mpz_t(), x.d, x.n);
		if (negative)
		{
			mpz_neg(host_x.get_mpz_t(), host_x.get_mpz_t());
		}
		mpz_fdiv_q_2exp(host_x.get_mpz_t(), host_x.get_mpz_t(), GMP_LIMB_BITS * limbs);
		mpz_abs(host_x.get_mpz_t(), host_x.get_mpz_t());
		From_mpz(temp2, host_x.get_mpz_t(), 2);
		out_limbs = temp2.n;
		return temp2.d;
	}
	if (negative)
	{
		gpu->Add(x.d, x.n, x.d, x.n, NULL, 0, 0, 1, true, NULL, 0);
	}
	uint64_t* result = x.d + limbs;
	out_limbs = x.n - limbs;
	if (negative)
	{
		gpu->Add(result, out_limbs, result, out_limbs, NULL, 0, 0, 1, false, NULL, 0);
	}
	return result;
}

//result = a modulo 2^L - 1, from 0 to 2^L - 2, as mod_2expm1 in Tree_backend.cpp.
void Gpu_root::Mod_2expm1(Gpu_number& result, const uint64_t* a, size_t a_limbs, mp_bitcnt_t L)
{
	const size_t L_limbs = (size_t)(L / 64) + 2;
	Zero(result, (a_limbs > L_limbs) ? a_limbs : L_limbs);
	CUDA_CHECK(cudaMemcpyAsync(result.d, a, a_limbs * sizeof(uint64_t), cudaMemcpyDeviceToDevice, 0));
	while (Bits(result.d, result.n) > L)
	{
		const size_t n = result.n;
		temp.Resize(n);
		temp2.Resize(n);
		gpu->Shift(temp.d, n, result.d, n, (long long)L, 0, false, 0);
		gpu->Shift(temp2.d, n, result.d, n, 0, L, false, 0);
		gpu->Add(result.d, n, temp2.d, n, temp.d, n, 0, 0, false, NULL, 0);
	}
	temp.Resize(result.n);
	gpu->Add(temp.d, temp.n, result.d, result.n, NULL, 0, 0, 1, false, NULL, 0);
	if (Bits(temp.d, temp.n) > L)
	{
		CUDA_CHECK(cudaMemsetAsync(result.d, 0, result.n * sizeof(uint64_t), 0));
	}
}

//r = a mod x, as Cpu_tree_backend::Barrett_reduce, with the inverse in level[0]. r must not be a.
void Gpu_root::Barrett_reduce(Gpu_number& r, const uint64_t* a, size_t a_limbs, mp_bitcnt_t a_bits, const uint64_t* x, size_t x_limbs, mp_bitcnt_t e, mp_bitcnt_t bits)
{
	if (a_bits < bits)
	{
		Zero(r, a_limbs + 1);
		CUDA_CHECK(cudaMemcpyAsync(r.d, a, a_limbs * sizeof(uint64_t), cudaMemcpyDeviceToDevice, 0));
		return;
	}

	const Gpu_number& V = level[0];
	const mp_bitcnt_t quotient_bits = a_bits - bits + 1;
	const mp_bitcnt_t cut = (bits > quotient_bits + GMP_LIMB_BITS) ? bits - quotient_bits - GMP_LIMB_BITS : 0;
	const mp_bitcnt_t v_shift = e - 2 * bits + cut, a_shift = bits - 1, q_shift = bits + 1 - cut;
	quotient.Resize((V.n > v_shift / 64 + 1) ? V.n - v_shift / 64 : 1);
	gpu->Shift(quotient.d, quotient.n, V.d, V.n, (long long)v_shift, 0, false, 0);
	temp.Resize((a_limbs > a_shift / 64 + 1) ? a_limbs - a_shift / 64 : 1);
	gpu->Shift(temp.d, temp.n, a, a_limbs, (long long)a_shift, 0, false, 0);
	Multiply(product, quotient.d, quotient.n, temp.d, temp.n);
	quotient.Resize((product.n > q_shift / 64 + 1) ? product.n - q_shift / 64 : 1);
	gpu->Shift(quotient.d, quotient.n, product.d, product.n, (long long)q_shift, 0, false, 0);

	//r = a - q x modulo 2^L - 1
	const mp_bitcnt_t L = Product_mod(product, quotient.d, quotient.n, x, x_limbs, bits + 3);
	Mod_2expm1(part, a, a_limbs, L);
	const size_t n = (part.n > product.n) ? part.n : product.n;
	r.Resize(n);
	gpu->Add(r.d, n, part.d, part.n, product.d, product.n, 0, 0, true, words + 6, 0);
	if (Carry())
	{
		gpu->Add(r.d, n, r.d, n, NULL, 0, (size_t)(L / 64), 1ULL << (L % 64), false, NULL, 0);
		gpu->Add(r.d, n, r.d, n, NULL, 0, 0, 1, true, NULL, 0);
	}
	for (int i = 0;; i++)
	{
		if (i == 64)
		{
			fprintf(stderr, "Gpu_root: the Barrett reduction did not converge\n");
			exit(3);
		}
		temp.Resize(n);
		gpu->Add(temp.d, n, r.d, n, x, x_limbs, 0, 0, true, words + 6, 0);
		if (Carry())
		{
			break; //r < x
		}
		r.Swap(temp);
	}
}

void Gpu_root::Compute(Gpu_number& fraction, const uint64_t* X, size_t x_limbs, mp_bitcnt_t bits, mp_bitcnt_t precision, mpz_srcptr prime_product)
{
	const size_t limbs = x_limbs;
	const mp_bitcnt_t normalize = limbs * GMP_LIMB_BITS - bits;

	//V, the inverse of the root shifted to fill limbs + 1 limbs, as in Cpu_tree_backend::Scaled_root.
	shifted.Resize(limbs + 1);
	gpu->Shift(shifted.d, shifted.n, X, x_limbs, -(long long)(normalize + GMP_LIMB_BITS), 0, false, 0);
	Newton_inverse(0, shifted.d, limbs + 1);
	if (limbs + 1 > INVERSE_EXACT_LIMBS)
	{
		gpu->Add(level[0].d, level[0].n, level[0].d, level[0].n, NULL, 0, 0, INVERSE_MARGIN, true, NULL, 0);
	}
	const mp_bitcnt_t e = 128 * limbs + GMP_LIMB_BITS - normalize;

	//The prime product is copied to the GPU once.
	const size_t p_limbs = mpz_size(prime_product);
	if (prime_product_host != prime_product || prime_product_limbs != p_limbs)
	{
		From_mpz(prime_product_gpu, prime_product, 1);
		prime_product_host = prime_product;
		prime_product_limbs = p_limbs;
	}
	const mp_bitcnt_t product_size = mpz_sizeinbase(prime_product, 2);

	if (product_size <= 2 * bits)
	{
		Barrett_reduce(remainder, prime_product_gpu.d, prime_product_gpu.n, product_size, X, x_limbs, e, bits);
	}
	else
	{
		//One part at a time from the top, as on the CPU.
		mp_bitcnt_t position = product_size;
		Zero(remainder, 1);
		while (position > 0)
		{
			const mp_bitcnt_t size = (position > bits - 1) ? bits - 1 : position;
			position -= size;
			//r = r 2^size + the next size bits of P
			const size_t n = remainder.n + (size_t)(size / 64) + 2;
			d.Resize(n);
			temp2.Resize(n);
			gpu->Shift(d.d, n, remainder.d, remainder.n, -(long long)size, 0, false, 0);
			gpu->Shift(temp2.d, n, prime_product_gpu.d, prime_product_gpu.n, (long long)position, size, false, 0);
			gpu->Add(d.d, n, d.d, n, temp2.d, n, 0, 0, false, NULL, 0);
			Barrett_reduce(remainder, d.d, d.n, Bits(d.d, d.n), X, x_limbs, e, bits);
		}
	}

	//floor(r 2^precision / X) = floor(r V / 2^(e - precision)), or one less.
	Multiply(product, remainder.d, remainder.n, level[0].d, level[0].n);
	fraction.Resize((size_t)((precision + 63) / 64) + 1);
	gpu->Shift(fraction.d, fraction.n, product.d, product.n, (long long)(e - precision), 0, false, 0);
}

class Gpu_tree_backend : public Cpu_tree_backend
{
public:
	Gpu_tree_backend(int levels, int max_bits_in_input);
	~Gpu_tree_backend();
	mp_bitcnt_t Build_product_tree();
	void Scaled_root(mpz_srcptr prime_product);
	void Remainder_tree();
	void Leaf_cofactors();

private:
	//The nodes of a depth: the product tree, and then the fractions of the remainder tree.
	struct Depth
	{
		Gpu_level level;
		size_t limb_capacity, size_capacity;
		uint32_t* bits; //of the product tree nodes
		size_t bits_capacity;
		uint32_t* transform; //of the nodes, kept from the product tree where the depth uses the NTT
		size_t transform_capacity;
		Ntt_plan plan;
		bool uses_ntt;
	};

	void Set_stride(Depth& d, size_t count, size_t stride);
	void Upload(int depth);
	void Download(int depth);
	void Download_bits(int depth);
	void Shifts(int depth);

	Gpu_ntt* gpu;
	Gpu_root* root;
	Gpu_number fraction; //of the root
	size_t leaf_limbs; //the most limbs of a leaf
	Gpu_level cofactors; //of the leaves
	size_t cofactor_limb_capacity, cofactor_size_capacity;
	std::vector<Depth> depth_data;
	uint32_t* parent_transform;
	size_t parent_transform_capacity;
	uint32_t* pair_transform;
	size_t pair_transform_capacity;
	Gpu_level scratch_level; //cyclic products, and exact truncated products
	size_t scratch_limb_capacity, scratch_size_capacity;
	uint32_t *shift, *precision;
	size_t shift_capacity, precision_capacity;
	std::vector<uint64_t> host_limbs;
	std::vector<uint32_t> host_sizes, host_words;
};

Gpu_tree_backend::Gpu_tree_backend(int levels, int max_bits_in_input) : Cpu_tree_backend(new Gpu_tree_products, levels, max_bits_in_input)
{
	gpu = &static_cast<Gpu_tree_products*>(products)->gpu;
	root = new Gpu_root(gpu);
	depth_data.resize(levels);
	for (int d = 0; d < levels; d++)
	{
		Depth& x = depth_data[d];
		x.level.limbs = NULL;
		x.level.size = NULL;
		x.level.stride = 0;
		x.level.count = (size_t)1 << d;
		x.limb_capacity = x.size_capacity = x.bits_capacity = x.transform_capacity = 0;
		x.bits = NULL;
		x.transform = NULL;
		x.uses_ntt = false;
		reserve(x.level.size, x.size_capacity, x.level.count);
		reserve(x.bits, x.bits_capacity, x.level.count);
	}
	parent_transform = pair_transform = NULL;
	parent_transform_capacity = pair_transform_capacity = 0;
	scratch_level.limbs = NULL;
	scratch_level.size = NULL;
	scratch_limb_capacity = scratch_size_capacity = 0;
	shift = precision = NULL;
	shift_capacity = precision_capacity = 0;
	leaf_limbs = 1;
	cofactors.limbs = NULL;
	cofactors.size = NULL;
	cofactors.stride = 0;
	cofactors.count = 0;
	cofactor_limb_capacity = cofactor_size_capacity = 0;
}

Gpu_tree_backend::~Gpu_tree_backend()
{
	delete root;
	for (Depth& x : depth_data)
	{
		cudaFree(x.level.limbs);
		cudaFree(x.level.size);
		cudaFree(x.bits);
		cudaFree(x.transform);
	}
	cudaFree(parent_transform);
	cudaFree(pair_transform);
	cudaFree(scratch_level.limbs);
	cudaFree(scratch_level.size);
	cudaFree(shift);
	cudaFree(precision);
	cudaFree(cofactors.limbs);
	cudaFree(cofactors.size);
}

void Gpu_tree_backend::Set_stride(Depth& d, size_t count, size_t stride)
{
	d.level.stride = stride;
	d.level.count = count;
	reserve(d.level.limbs, d.limb_capacity, count * stride);
}

//Copy the nodes of a depth from the CPU to the GPU.
void Gpu_tree_backend::Upload(int depth)
{
	Depth& d = depth_data[depth];
	const int first = (1 << depth) - 1;
	const size_t count = d.level.count, stride = d.level.stride;
	host_limbs.assign(count * stride, 0);
	host_sizes.resize(count);
	for (size_t r = 0; r < count; r++)
	{
		mpz_srcptr a = node[first + r].get_mpz_t();
		size_t n = mpz_size(a);
		if (n > stride)
		{
			fprintf(stderr, "Gpu_tree_backend: a node of %zu limbs does not fit in %zu\n", n, stride);
			exit(3);
		}
		memcpy(host_limbs.data() + r * stride, a->_mp_d, n * sizeof(uint64_t));
		host_sizes[r] = (uint32_t)n;
	}
	CUDA_CHECK(cudaMemcpy(d.level.limbs, host_limbs.data(), host_limbs.size() * sizeof(uint64_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(d.level.size, host_sizes.data(), count * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

//Copy the numbers of a depth from the GPU to the CPU.
void Gpu_tree_backend::Download(int depth)
{
	Depth& d = depth_data[depth];
	const int first = (1 << depth) - 1;
	const size_t count = d.level.count, stride = d.level.stride;
	host_limbs.resize(count * stride);
	host_sizes.resize(count);
	CUDA_CHECK(cudaMemcpy(host_limbs.data(), d.level.limbs, host_limbs.size() * sizeof(uint64_t), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(host_sizes.data(), d.level.size, count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
	for (size_t r = 0; r < count; r++)
	{
		mpz_ptr a = node[first + r].get_mpz_t();
		size_t n = host_sizes[r];
		if (a->_mp_alloc < (int)n)
		{
			mpz_realloc2(a, n * GMP_LIMB_BITS);
		}
		memcpy(a->_mp_d, host_limbs.data() + r * stride, n * sizeof(uint64_t));
		a->_mp_size = (int)n;
	}
}

void Gpu_tree_backend::Download_bits(int depth)
{
	Depth& d = depth_data[depth];
	const int first = (1 << depth) - 1;
	host_words.resize(d.level.count);
	gpu->Bits(d.level, d.bits, 0);
	CUDA_CHECK(cudaMemcpy(host_words.data(), d.bits, d.level.count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
	for (size_t r = 0; r < d.level.count; r++)
	{
		node_bits[first + r] = host_words[r];
	}
}

mp_bitcnt_t Gpu_tree_backend::Build_product_tree()
{
	//The leaves, with room for their fractions too.
	const int leaf_depth = levels - 1;
	size_t max_leaf = 1;
	for (int i = first_leaf; i < first_leaf + num_leaves; i++)
	{
		node_bits[i] = mpz_sizeinbase(node[i].get_mpz_t(), 2);
		size_t n = mpz_size(node[i].get_mpz_t());
		max_leaf = (n > max_leaf) ? n : max_leaf;
	}
	Set_stride(depth_data[leaf_depth], num_leaves, max_leaf + 1);
	Upload(leaf_depth);
	gpu->Bits(depth_data[leaf_depth].level, depth_data[leaf_depth].bits, 0);
	leaf_limbs = max_leaf;

	for (int depth = levels - 2; depth >= 0; depth--)
	{
		Depth& children = depth_data[depth + 1];
		Depth& parents = depth_data[depth];
		const size_t count = parents.level.count;
		children.uses_ntt = tree_level_plan(children.plan, &node_bits[(1 << (depth + 1)) - 1], children.level.count, guard);
		if (children.uses_ntt)
		{
			const Ntt_plan& p = children.plan;
			const size_t T = p.Length();
			//The transforms of the children are kept for the remainder tree, except those of the leaves.
			uint32_t* X;
			if (depth + 1 <= levels - 2)
			{
				reserve(children.transform, children.transform_capacity, (size_t)p.primes * 2 * count * T);
				X = children.transform;
			}
			else
			{
				reserve(pair_transform, pair_transform_capacity, (size_t)p.primes * 2 * count * T);
				X = pair_transform;
			}
			reserve(parent_transform, parent_transform_capacity, (size_t)p.primes * count * T);
			Set_stride(parents, count, Gpu_ntt::Combine_limbs(p));
			gpu->Forward(p, children.level, X, 0);
			gpu->Pointwise_adjacent(p, X, parent_transform, count, 0);
			gpu->Combine(p, parent_transform, count, false, parents.level, 0);
		}
		else
		{
			Set_stride(parents, count, 2 * children.level.stride);
			gpu->School_products(children.level, parents.level, 0);
		}
		Download_bits(depth);
	}
	return node_bits[0];
}

//The root stays on the GPU: the fraction of the root replaces the product there.
void Gpu_tree_backend::Scaled_root(mpz_srcptr prime_product)
{
	Depth& top = depth_data[0];
	const size_t limbs = (size_t)((node_bits[0] + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS);
	root->Compute(fraction, top.level.limbs, limbs, node_bits[0], node_bits[0] + guard, prime_product);
	if (fraction.n > top.level.stride)
	{
		Set_stride(top, 1, fraction.n);
	}
	CUDA_CHECK(cudaMemcpy(top.level.limbs, fraction.d, fraction.n * sizeof(uint64_t), cudaMemcpyDeviceToDevice));
	gpu->Size(top.level.limbs, fraction.n, top.level.size, NULL, 0);
	if (levels < 3)
	{
		Download(0); //The leaves are below the root.
	}
}

//The shifts and precisions of the middle products of the nodes of a depth: precision = parent bits + guard, and
//shift = parent bits - node bits.
void Gpu_tree_backend::Shifts(int depth)
{
	const int first = (1 << depth) - 1;
	const size_t count = (size_t)1 << depth;
	host_words.resize(2 * count);
	for (size_t r = 0; r < count; r++)
	{
		const int i = first + (int)r, parent = (i - 1) / 2;
		host_words[r] = (uint32_t)(node_bits[parent] - node_bits[i]);
		host_words[count + r] = (uint32_t)(node_bits[parent] + guard);
	}
	reserve(shift, shift_capacity, count);
	reserve(precision, precision_capacity, count);
	CUDA_CHECK(cudaMemcpy(shift, host_words.data(), count * sizeof(uint32_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(precision, host_words.data() + count, count * sizeof(uint32_t), cudaMemcpyHostToDevice));
}

void Gpu_tree_backend::Remainder_tree()
{
	if (levels < 3)
	{
		return; //The root is the parent of the leaves, and is on the CPU.
	}
	for (int depth = 1; depth <= levels - 2; depth++)
	{
		Depth& nodes = depth_data[depth];
		Depth& parents = depth_data[depth - 1];
		const size_t count = nodes.level.count;
		Shifts(depth);
		if (nodes.uses_ntt)
		{
			const Ntt_plan& p = nodes.plan;
			const size_t T = p.Length();
			reserve(parent_transform, parent_transform_capacity, (size_t)p.primes * (count / 2) * T);
			reserve(pair_transform, pair_transform_capacity, (size_t)p.primes * count * T);
			const size_t combine_limbs = Gpu_ntt::Combine_limbs(p);
			reserve(scratch_level.limbs, scratch_limb_capacity, count * combine_limbs);
			reserve(scratch_level.size, scratch_size_capacity, count);
			scratch_level.stride = combine_limbs;
			scratch_level.count = count;
			gpu->Forward(p, parents.level, parent_transform, 0);
			gpu->Pointwise_swap(p, parent_transform, nodes.transform, pair_transform, count / 2, 0);
			gpu->Combine(p, pair_transform, count, true, scratch_level, 0);
			gpu->Extract_middle(scratch_level, shift, precision, parents.level.size, 2, nodes.level, 0);
		}
		else
		{
			//The exact truncated products need the nodes themselves, so they are written to the scratch level, which then
			//takes the place of the nodes.
			reserve(scratch_level.limbs, scratch_limb_capacity, count * nodes.level.stride);
			reserve(scratch_level.size, scratch_size_capacity, count);
			scratch_level.stride = nodes.level.stride;
			scratch_level.count = count;
			gpu->School_middle(parents.level, nodes.level, shift, precision, scratch_level, 0);
			std::swap(nodes.level.limbs, scratch_level.limbs);
			std::swap(nodes.level.size, scratch_level.size);
			std::swap(nodes.limb_capacity, scratch_limb_capacity);
			std::swap(nodes.size_capacity, scratch_size_capacity);
		}
	}
}

void Gpu_tree_backend::Leaf_cofactors()
{
	Depth& leaves = depth_data[levels - 1];
	Depth& parents = depth_data[levels - 2];
	const int K = Gpu_ntt::Leaf_kernel_limbs((int)leaf_limbs);
	if (K != 0)
	{
		reserve(cofactors.limbs, cofactor_limb_capacity, (size_t)num_leaves * K);
		reserve(cofactors.size, cofactor_size_capacity, (size_t)num_leaves);
		cofactors.stride = K;
		cofactors.count = num_leaves;
		if (gpu->Leaf_cofactors(parents.level, parents.bits, leaves.level, leaves.bits, (int)leaf_limbs, (uint32_t)guard, MAX_SQUARINGS_IN_CHECK, cofactors, 0))
		{
			host_limbs.resize((size_t)num_leaves * K);
			host_sizes.resize(num_leaves);
			CUDA_CHECK(cudaMemcpy(host_limbs.data(), cofactors.limbs, host_limbs.size() * sizeof(uint64_t), cudaMemcpyDeviceToHost));
			CUDA_CHECK(cudaMemcpy(host_sizes.data(), cofactors.size, num_leaves * sizeof(uint32_t), cudaMemcpyDeviceToHost));
			parallel_for(num_leaves, 4096, [&](int r)
			{
				mpz_ptr a = cofactor[r].get_mpz_t();
				const int n = (int)host_sizes[r];
				if (a->_mp_alloc < n)
				{
					mpz_realloc2(a, n * GMP_LIMB_BITS);
				}
				memcpy(a->_mp_d, host_limbs.data() + (size_t)r * K, n * sizeof(uint64_t));
				a->_mp_size = n;
			});
			return;
		}
	}
	//Large leaves: on the CPU, from the fractions of their parents.
	Download(levels - 2);
	Cpu_tree_backend::Leaf_cofactors();
}

//Before the first use of the GPU: the threads that wait for it sleep instead of spinning, which leaves the CPU to the threads
//that check the batches (the waits get a little longer).
static void gpu_start()
{
	static bool started = false;
	if (!started)
	{
		started = true;
		cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);
	}
}

//The product of n factors on the GPU, for the prime products: a product tree whose leaves are products of up to
//PRODUCT_LEAF_FACTORS factors, computed on the CPU, with one batched step per depth as in Build_product_tree: schoolbook
//products for small nodes and NTT products above. A depth beyond the longest transform is finished on the CPU.
#define PRODUCT_LEAF_FACTORS 16
void gpu_product_list(const uint32_t* factors, size_t n, mpz_class& result)
{
	gpu_start();
	const size_t leaves = (n + PRODUCT_LEAF_FACTORS - 1) / PRODUCT_LEAF_FACTORS;
	size_t count = 2;
	while (count < leaves)
	{
		count *= 2;
	}
	const size_t leaf_stride = PRODUCT_LEAF_FACTORS / 2 + 1;
	std::vector<uint64_t> host_limbs(count * leaf_stride, 0);
	std::vector<uint32_t> host_sizes(count, 1);
	for (size_t r = 0; r < count; r++)
	{
		mp_limb_t* a = (mp_limb_t*)host_limbs.data() + r * leaf_stride;
		mp_size_t size = 1;
		a[0] = 1;
		for (size_t k = r * PRODUCT_LEAF_FACTORS; k < n && k < (r + 1) * PRODUCT_LEAF_FACTORS; k++)
		{
			mp_limb_t carry = mpn_mul_1(a, a, size, factors[k]);
			if (carry != 0)
			{
				a[size++] = carry;
			}
		}
		host_sizes[r] = (uint32_t)size;
	}

	Gpu_ntt gpu;
	Gpu_level level[2];
	size_t limb_capacity[2] = { 0, 0 }, size_capacity[2] = { 0, 0 }, bits_capacity = 0, x_capacity = 0, y_capacity = 0;
	uint32_t *bits = NULL, *X = NULL, *Y = NULL;
	for (int k = 0; k < 2; k++)
	{
		level[k].limbs = NULL;
		level[k].size = NULL;
	}
	reserve(level[0].limbs, limb_capacity[0], host_limbs.size());
	reserve(level[0].size, size_capacity[0], count);
	level[0].stride = leaf_stride;
	level[0].count = count;
	CUDA_CHECK(cudaMemcpy(level[0].limbs, host_limbs.data(), host_limbs.size() * sizeof(uint64_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(level[0].size, host_sizes.data(), count * sizeof(uint32_t), cudaMemcpyHostToDevice));

	int cur = 0;
	while (level[cur].count > 1)
	{
		Gpu_level& children = level[cur];
		Gpu_level& parents = level[1 - cur];
		const size_t pairs = children.count / 2;
		reserve(bits, bits_capacity, children.count);
		gpu.Bits(children, bits, 0);
		std::vector<uint32_t> host_bits(children.count);
		CUDA_CHECK(cudaMemcpy(host_bits.data(), bits, children.count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
		const uint32_t max_bits = *std::max_element(host_bits.begin(), host_bits.end());
		reserve(parents.size, size_capacity[1 - cur], pairs);
		parents.count = pairs;
		Ntt_plan plan;
		if (max_bits < TREE_NTT_MIN_CHILD_BITS)
		{
			parents.stride = 2 * children.stride;
			reserve(parents.limbs, limb_capacity[1 - cur], pairs * parents.stride);
			gpu.School_products(children, parents, 0);
		}
		else if (ntt_plan_product(plan, max_bits, max_bits))
		{
			const size_t T = plan.Length();
			reserve(X, x_capacity, (size_t)plan.primes * children.count * T);
			reserve(Y, y_capacity, (size_t)plan.primes * pairs * T);
			parents.stride = Gpu_ntt::Combine_limbs(plan);
			reserve(parents.limbs, limb_capacity[1 - cur], pairs * parents.stride);
			gpu.Forward(plan, children, X, 0);
			gpu.Pointwise_adjacent(plan, X, Y, pairs, 0);
			gpu.Combine(plan, Y, pairs, false, parents, 0);
		}
		else
		{
			break; //Beyond the longest transform.
		}
		cur = 1 - cur;
	}

	//The numbers of the last depth, multiplied on the CPU if there is more than one.
	Gpu_level& top = level[cur];
	host_limbs.resize(top.count * top.stride);
	host_sizes.resize(top.count);
	CUDA_CHECK(cudaMemcpy(host_limbs.data(), top.limbs, host_limbs.size() * sizeof(uint64_t), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(host_sizes.data(), top.size, top.count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
	std::vector<mpz_class> part(top.count);
	for (size_t r = 0; r < top.count; r++)
	{
		mpz_import(part[r].get_mpz_t(), host_sizes[r], -1, sizeof(uint64_t), 0, 0, host_limbs.data() + r * top.stride);
	}
	for (size_t width = 1; width < top.count; width *= 2)
	{
		for (size_t r = 0; r + width < top.count; r += 2 * width)
		{
			part[r] *= part[r + width];
		}
	}
	result = part[0];

	for (int k = 0; k < 2; k++)
	{
		cudaFree(level[k].limbs);
		cudaFree(level[k].size);
	}
	cudaFree(bits);
	cudaFree(X);
	cudaFree(Y);
}

Tree_backend* create_gpu_tree_backend(int levels, int max_bits_in_input)
{
	gpu_start();
	return new Gpu_tree_backend(levels, max_bits_in_input);
}
