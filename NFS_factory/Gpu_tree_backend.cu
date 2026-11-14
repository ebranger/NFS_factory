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
#include <mutex>
#include <condition_variable>
#include <vector>

static void gpu_start();

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
		gpu.Combine_products(plan, PRODUCTS_ADJACENT, NULL, X, Y, 1, cyclic, out, 0);
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
//With some room to spare, since the sizes vary a little from batch to batch, and cudaFree waits for the whole GPU.
template <class T> static void reserve(T*& p, size_t& capacity, size_t count)
{
	if (count > capacity)
	{
		cudaFree(p);
		capacity = count + count / 8;
		CUDA_CHECK(cudaMalloc(&p, capacity * sizeof(T)));
	}
}

//Page-locked host memory for the large copies to and from the GPU, which then run at full speed. It grows as needed and is kept;
//resize does not keep the contents.
template <class T> class Pinned
{
public:
	Pinned() : p(NULL), capacity(0), count(0) {}
	~Pinned() { cudaFreeHost(p); }
	Pinned(const Pinned&) = delete;
	Pinned& operator=(const Pinned&) = delete;
	void resize(size_t n)
	{
		if (n > capacity)
		{
			cudaFreeHost(p);
			capacity = n + n / 8;
			CUDA_CHECK(cudaMallocHost(&p, capacity * sizeof(T)));
		}
		count = n;
	}
	T* data() { return p; }
	size_t size() const { return count; }
	T& operator[](size_t i) { return p[i]; }

private:
	T* p;
	size_t capacity, count;
};

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
#define ROOT_CACHE_ALIGN 18 //the grain in bits of the sizes for which the transform of the top of the prime product is kept
#define PRIME_PRODUCT_KEY (1LL << 62) //the key of the transform of the prime product itself, kept by Gpu_root::Multiply

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
	void Barrett_reduce(Gpu_number& r, const uint64_t* a, size_t a_limbs, mp_bitcnt_t a_bits, const uint64_t* x, size_t x_limbs, mp_bitcnt_t e, mp_bitcnt_t bits,
		bool constant_a = false);
	void Mod_2expm1(Gpu_number& result, const uint64_t* a, size_t a_limbs, mp_bitcnt_t L);
	void Multiply(Gpu_number& out, const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs, long long b_key = -1);
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
	//The transform of the b of the last Multiply with a b_key, which the next one with the same key and plan reuses.
	uint32_t* cached;
	size_t cached_words;
	long long cached_key;
	Ntt_plan cached_plan;
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
	cached = NULL;
	cached_words = 0;
	cached_key = -1;
	prime_product_host = NULL;
	prime_product_limbs = 0;
}

Gpu_root::~Gpu_root()
{
	cudaFree(words);
	cudaFree(cached);
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
		transform_words = words_needed + words_needed / 8;
		CUDA_CHECK(cudaMalloc(&X, transform_words * sizeof(uint32_t)));
		CUDA_CHECK(cudaMalloc(&Y, transform_words * sizeof(uint32_t)));
	}
	Gpu_level la = { const_cast<uint64_t*>(a), host_words[0], words, 1 };
	Gpu_level lb = { const_cast<uint64_t*>(b), host_words[2], words + 2, 1 };
	gpu->Forward(plan, la, X, 0);
	gpu->Forward(plan, lb, Y, 0);
	return true;
}

//out = a * b. out must not be a or b. With b_key (0 or more), b is the same number whenever the key is, and its transform is
//kept for the next product with the key and the same plan.
void Gpu_root::Multiply(Gpu_number& out, const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs, long long b_key)
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
	if (b_key >= 0)
	{
		//A plan for sizes rounded up, so that it is the same for products of about the same size, and b's transform can be kept.
		const mp_bitcnt_t grain = (mp_bitcnt_t)1 << ROOT_CACHE_ALIGN;
		Ntt_plan rounded;
		if (ntt_plan_product(rounded, (host_words[1] + grain - 1) / grain * grain, (host_words[3] + grain - 1) / grain * grain))
		{
			plan = rounded;
		}
	}
	out.Resize(Gpu_ntt::Combine_limbs(plan));
	Gpu_level o = { out.d, out.n, words + 4, 1 };
	if (b_key >= 0)
	{
		const bool same = (b_key == cached_key && plan.primes == cached_plan.primes && plan.log_length == cached_plan.log_length &&
			plan.coefficient_bits == cached_plan.coefficient_bits);
		const size_t words_needed = (size_t)plan.primes * plan.Length();
		if (!same)
		{
			if (words_needed > cached_words)
			{
				cudaFree(cached);
				cached_words = words_needed + words_needed / 8;
				CUDA_CHECK(cudaMalloc(&cached, cached_words * sizeof(uint32_t)));
			}
			Gpu_level lb = { const_cast<uint64_t*>(b), host_words[2], words + 2, 1 };
			gpu->Forward(plan, lb, cached, 0);
			cached_key = b_key;
			cached_plan = plan;
		}
		if (words_needed > transform_words)
		{
			cudaFree(X);
			cudaFree(Y);
			transform_words = words_needed + words_needed / 8;
			CUDA_CHECK(cudaMalloc(&X, transform_words * sizeof(uint32_t)));
			CUDA_CHECK(cudaMalloc(&Y, transform_words * sizeof(uint32_t)));
		}
		Gpu_level la = { const_cast<uint64_t*>(a), host_words[0], words, 1 };
		gpu->Forward(plan, la, X, 0);
		gpu->Combine_products(plan, PRODUCTS_PLAIN, X, cached, X, 1, false, o, 0);
		return;
	}
	Transform_pair(plan, a, b);
	gpu->Combine_products(plan, PRODUCTS_PLAIN, X, Y, X, 1, false, o, 0);
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
	gpu->Combine_products(plan, PRODUCTS_PLAIN, X, Y, X, 1, true, o, 0);
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
//With constant_a (a is the prime product), the top part of a that the quotient is estimated from starts at a multiple of
//2^ROOT_CACHE_ALIGN bits, so that it is the same number for roots of about the same size, and Multiply keeps its transform. The
//estimate then uses a few more bits of a, and is still not above the quotient.
void Gpu_root::Barrett_reduce(Gpu_number& r, const uint64_t* a, size_t a_limbs, mp_bitcnt_t a_bits, const uint64_t* x, size_t x_limbs, mp_bitcnt_t e, mp_bitcnt_t bits,
	bool constant_a)
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
	const mp_bitcnt_t extra = constant_a ? (bits - 1) & (((mp_bitcnt_t)1 << ROOT_CACHE_ALIGN) - 1) : 0;
	const mp_bitcnt_t v_shift = e - 2 * bits + cut, a_shift = bits - 1 - extra, q_shift = bits + 1 - cut + extra;
	quotient.Resize((V.n > v_shift / 64 + 1) ? V.n - v_shift / 64 : 1);
	gpu->Shift(quotient.d, quotient.n, V.d, V.n, (long long)v_shift, 0, false, 0);
	temp.Resize((a_limbs > a_shift / 64 + 1) ? a_limbs - a_shift / 64 : 1);
	gpu->Shift(temp.d, temp.n, a, a_limbs, (long long)a_shift, 0, false, 0);
	Multiply(product, quotient.d, quotient.n, temp.d, temp.n, constant_a ? (long long)a_shift : -1);
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
		cached_key = -1;
	}
	const mp_bitcnt_t product_size = mpz_sizeinbase(prime_product, 2);

	if (product_size < bits)
	{
		//The prime product is its own remainder, and is used as it is below.
	}
	else if (product_size <= 2 * bits)
	{
		Barrett_reduce(remainder, prime_product_gpu.d, prime_product_gpu.n, product_size, X, x_limbs, e, bits, true);
	}
	else if (product_size / CHUNKED_REDUCTION_MAX_PARTS >= bits)
	{
		//A root much smaller than the prime product (a batch of few relations, at the end of a file or of a second stage): the
		//reduction below would take a part of the size of the root at a time, each with tens of kernel launches. The remainder is
		//the same number however it is computed, so divide on the CPU.
		mpz_class x, r;
		To_mpz(x.get_mpz_t(), X, x_limbs);
		mpz_tdiv_r(r.get_mpz_t(), prime_product, x.get_mpz_t());
		From_mpz(remainder, r.get_mpz_t(), 1);
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

	//floor(r 2^precision / X) = floor(r V / 2^(e - precision)), or one less. With P less than X, r is P, and its transform is kept
	//for the next batches (the key is not a shift, which the Barrett reduction uses).
	if (product_size < bits)
	{
		Multiply(product, level[0].d, level[0].n, prime_product_gpu.d, prime_product_gpu.n, PRIME_PRODUCT_KEY);
	}
	else
	{
		Multiply(product, remainder.d, remainder.n, level[0].d, level[0].n);
	}
	fraction.Resize((size_t)((precision + 63) / 64) + 1);
	gpu->Shift(fraction.d, fraction.n, product.d, product.n, (long long)(e - precision), 0, false, 0);
}

static unsigned int grid_size_leaves(size_t count)
{
	const size_t blocks = (count + 127) / 128;
	return (unsigned int)((blocks < 65536) ? ((blocks > 0) ? blocks : 1) : 65536);
}

//The checks of the cofactors on the GPU, one thread per leaf: Batch_smooth::Classify_cofactor for cofactors of up to two limbs,
//with the same primality tests and the same Pollard-Brent rho as Smooth_factor.cpp, so that the results are the same. Cofactors
//that need more (the rare factoring of a one large prime cofactor, and cofactors of more than 128 bits in a range that is checked)
//are left to the CPU.
namespace
{

#define RHO_MAX_ITERATIONS (1 << 20) //as in Smooth_factor.cpp
#define RHO_GCD_INTERVAL 64

__device__ __forceinline__ uint64_t montgomery_inverse64(uint64_t n)
{
	uint64_t inverse = n; //correct to 3 bits, since n * n = 1 mod 8 for odd n.
	for (int i = 0; i < 5; i++)
	{
		inverse *= 2 - n * inverse;
	}
	return 0 - inverse;
}

//a b / 2^64 mod n, for a, b < n and odd n.
__device__ __forceinline__ uint64_t montgomery_multiply64(uint64_t a, uint64_t b, uint64_t n, uint64_t ninv)
{
	const uint64_t lo = a * b, hi = __umul64hi(a, b);
	const uint64_t m = lo * ninv;
	const uint64_t mhi = __umul64hi(m, n);
	//(lo + m n) is 0 mod 2^64, with a carry unless lo is 0.
	uint64_t t = hi + mhi;
	uint32_t carry = (t < hi);
	const uint64_t t2 = t + (lo != 0);
	carry |= (t2 < t);
	t = t2;
	if (carry || t >= n)
	{
		t -= n;
	}
	return t;
}

__device__ __forceinline__ uint64_t rho_step64(uint64_t y, uint64_t c, uint64_t n, uint64_t ninv)
{
	const uint64_t r = montgomery_multiply64(y, y, n, ninv);
	uint64_t sum = r + c;
	if (sum < r || sum >= n)
	{
		sum -= n;
	}
	return sum;
}

//The binary algorithm, since 64-bit division is slow on the GPU. The same gcd as in Smooth_factor.cpp.
__device__ uint64_t gcd64(uint64_t a, uint64_t b)
{
	if (a == 0 || b == 0)
	{
		return a | b;
	}
	const int shift = __ffsll((long long)(a | b)) - 1;
	a >>= __ffsll((long long)a) - 1;
	while (b != 0)
	{
		b >>= __ffsll((long long)b) - 1;
		if (a > b)
		{
			const uint64_t t = a;
			a = b;
			b = t;
		}
		b -= a;
	}
	return a << shift;
}

//Pollard-Brent rho for an odd composite n < 2^64, as rho64 in Smooth_factor.cpp. Returns 0 if no factor was found.
__device__ uint64_t rho64(uint64_t n)
{
	const uint64_t ninv = montgomery_inverse64(n);
	for (uint64_t c = 1; c <= 5; c++)
	{
		uint64_t y = 2, x = y, saved_y = y, q = 1, g = 1;
		uint64_t r = 1;
		while (g == 1 && r <= RHO_MAX_ITERATIONS)
		{
			x = y;
			for (uint64_t i = 0; i < r; i++)
			{
				y = rho_step64(y, c, n, ninv);
			}
			for (uint64_t k = 0; k < r && g == 1; k += RHO_GCD_INTERVAL)
			{
				saved_y = y;
				const uint64_t steps = (r - k < RHO_GCD_INTERVAL) ? r - k : RHO_GCD_INTERVAL;
				for (uint64_t i = 0; i < steps; i++)
				{
					y = rho_step64(y, c, n, ninv);
					q = montgomery_multiply64(q, (x > y) ? x - y : y - x, n, ninv);
				}
				g = gcd64(q, n);
			}
			r *= 2;
		}
		if (g == n || g == 0)
		{
			//Several factors were found at once, go back and take one step at a time.
			do
			{
				saved_y = rho_step64(saved_y, c, n, ninv);
				g = gcd64((x > saved_y) ? x - saved_y : saved_y - x, n);
			} while (g == 1);
		}
		if (g != 1 && g != n)
		{
			return g;
		}
	}
	return 0;
}

//Deterministic Miller-Rabin for n < 2^32, with the bases 2, 7 and 61.
__device__ bool is_prime32(uint64_t n)
{
	if (n < 2)
	{
		return false;
	}
	const uint64_t bases[3] = { 2, 7, 61 };
	for (int b = 0; b < 3; b++)
	{
		if (n == bases[b])
		{
			return true;
		}
		if (n % bases[b] == 0)
		{
			return false;
		}
	}
	uint64_t d = n - 1;
	int s = 0;
	while ((d & 1) == 0)
	{
		d >>= 1;
		s++;
	}
	for (int b = 0; b < 3; b++)
	{
		//x = base^d mod n
		uint64_t x = 1, base = bases[b] % n, e = d;
		while (e > 0)
		{
			if (e & 1)
			{
				x = x * base % n;
			}
			base = base * base % n;
			e >>= 1;
		}
		if (x == 1 || x == n - 1)
		{
			continue;
		}
		bool composite = true;
		for (int i = 1; i < s && composite; i++)
		{
			x = x * x % n;
			composite = (x != n - 1);
		}
		if (composite)
		{
			return false;
		}
	}
	return true;
}

//Strong probable prime to base 2, for odd n, 2^32 <= n < 2^64.
__device__ bool sprp2_64(uint64_t n)
{
	const uint64_t ninv = montgomery_inverse64(n);
	const uint64_t one = (0 - n) % n;
	const uint64_t minus_one = n - one;
	const int s = __ffsll((long long)(n - 1)) - 1;
	const uint64_t d = (n - 1) >> s;
	const int top = 63 - __clzll((long long)d);
	uint64_t x = one;
	for (int i = top; i >= 0; i--)
	{
		x = montgomery_multiply64(x, x, n, ninv);
		if ((d >> i) & 1)
		{
			const uint64_t doubled = x + x;
			x = (doubled < x || doubled >= n) ? doubled - n : doubled;
		}
	}
	if (x == one || x == minus_one)
	{
		return true;
	}
	for (int i = 1; i < s; i++)
	{
		x = montgomery_multiply64(x, x, n, ninv);
		if (x == minus_one)
		{
			return true;
		}
	}
	return false;
}

//x = r mod n for r < 2n (r of 129 bits, with carry the top bit).
__device__ __forceinline__ void reduce_once128(uint64_t& lo, uint64_t& hi, uint64_t carry, uint64_t nlo, uint64_t nhi)
{
	const uint64_t slo = lo - nlo;
	const uint64_t b1 = (lo < nlo);
	const uint64_t shi = hi - nhi - b1;
	const uint64_t borrow = (hi < nhi) | ((hi == nhi) & b1);
	if (carry || !borrow)
	{
		lo = slo;
		hi = shi;
	}
}

__device__ __forceinline__ void add_carry(uint64_t& sum, uint64_t a, uint64_t& carry)
{
	const uint64_t s = sum + a;
	const uint64_t c1 = (s < a);
	const uint64_t s2 = s + carry;
	carry = c1 + (s2 < s);
	sum = s2;
}

//a^2 / 2^128 mod n for a < n and odd n of two limbs, as montgomery_square in Smooth_factor.cpp.
__device__ void montgomery_square128(uint64_t& alo, uint64_t& ahi, uint64_t nlo, uint64_t nhi, uint64_t ninv)
{
	//t = a^2 = alo^2 + 2 alo ahi 2^64 + ahi^2 2^128, less than 2^256.
	uint64_t t0 = alo * alo, t1 = __umul64hi(alo, alo);
	uint64_t t2 = ahi * ahi, t3 = __umul64hi(ahi, ahi);
	uint64_t lo = alo * ahi, hi = __umul64hi(alo, ahi);
	const uint64_t top = hi >> 63;
	hi = (hi << 1) | (lo >> 63);
	lo <<= 1;
	uint64_t c = 0;
	add_carry(t1, lo, c);
	add_carry(t2, hi, c);
	add_carry(t3, top, c);
	//Two reduction steps, each clearing the lowest limb. c4 is the carry out of t3.
	uint64_t c4 = 0;
	uint64_t m = t0 * ninv;
	lo = m * nlo;
	hi = __umul64hi(m, nlo);
	c = 0;
	add_carry(t0, lo, c);
	uint64_t h = hi;
	lo = m * nhi;
	hi = __umul64hi(m, nhi);
	add_carry(t1, h, c);
	hi += c;
	c = 0;
	add_carry(t1, lo, c);
	add_carry(t2, hi, c);
	add_carry(t3, 0, c);
	c4 += c;
	m = t1 * ninv;
	lo = m * nlo;
	hi = __umul64hi(m, nlo);
	c = 0;
	add_carry(t1, lo, c);
	h = hi;
	lo = m * nhi;
	hi = __umul64hi(m, nhi);
	add_carry(t2, h, c);
	hi += c;
	c = 0;
	add_carry(t2, lo, c);
	add_carry(t3, hi, c);
	c4 += c;
	alo = t2;
	ahi = t3;
	reduce_once128(alo, ahi, c4, nlo, nhi);
}

__device__ __forceinline__ void double_mod128(uint64_t& lo, uint64_t& hi, uint64_t nlo, uint64_t nhi)
{
	const uint64_t carry = hi >> 63;
	hi = (hi << 1) | (lo >> 63);
	lo <<= 1;
	reduce_once128(lo, hi, carry, nlo, nhi);
}

//Strong probable prime to base 2, for odd n, 2^64 <= n < 2^128.
__device__ bool sprp2_128(uint64_t nlo, uint64_t nhi)
{
	const uint64_t ninv = montgomery_inverse64(nlo);
	//one = 2^128 mod n, by doubling 2^(bits - 1) < n up to 2^128.
	const int top = 63 - __clzll((long long)nhi);
	uint64_t one_lo = 0, one_hi = 1ULL << top;
	for (int i = top + 64; i < 128; i++)
	{
		double_mod128(one_lo, one_hi, nlo, nhi);
	}
	const uint64_t minus_lo = nlo - one_lo;
	const uint64_t minus_hi = nhi - one_hi - (nlo < one_lo);

	//n - 1 = d 2^s
	uint64_t dlo = nlo - 1, dhi = nhi;
	int s;
	if (dlo != 0)
	{
		s = __ffsll((long long)dlo) - 1;
		dlo = (dlo >> s) | ((s != 0) ? dhi << (64 - s) : 0);
		dhi >>= s;
	}
	else
	{
		s = 64 + __ffsll((long long)dhi) - 1;
		dlo = dhi >> (s - 64);
		dhi = 0;
	}
	const int d_top = (dhi != 0) ? 127 - __clzll((long long)dhi) : 63 - __clzll((long long)dlo);

	uint64_t xlo = one_lo, xhi = one_hi;
	for (int i = d_top; i >= 0; i--)
	{
		montgomery_square128(xlo, xhi, nlo, nhi, ninv);
		const uint64_t bit = (i >= 64) ? (dhi >> (i - 64)) & 1 : (dlo >> i) & 1;
		if (bit)
		{
			double_mod128(xlo, xhi, nlo, nhi);
		}
	}
	if ((xlo == one_lo && xhi == one_hi) || (xlo == minus_lo && xhi == minus_hi))
	{
		return true;
	}
	for (int i = 1; i < s; i++)
	{
		montgomery_square128(xlo, xhi, nlo, nhi, ninv);
		if (xlo == minus_lo && xhi == minus_hi)
		{
			return true;
		}
	}
	return false;
}

//is_probable_prime of Smooth_factor.cpp, for numbers of up to two limbs.
__device__ bool is_probable_prime128(uint64_t lo, uint64_t hi)
{
	if ((lo & 1) == 0)
	{
		return hi == 0 && lo == 2;
	}
	if (hi != 0)
	{
		return sprp2_128(lo, hi);
	}
	if (lo < (1ULL << 32))
	{
		return is_prime32(lo);
	}
	return sprp2_64(lo);
}

__device__ __forceinline__ int bits64(uint64_t x)
{
	return (x == 0) ? 0 : 64 - __clzll((long long)x);
}

//Check_two_primes of Batch_smooth.
__device__ unsigned char check_two_primes(uint64_t p, uint64_t q, const Leaf_check& check, uint64_t* lp)
{
	if (bits64(p) > check.one_lp_bits || bits64(q) > check.one_lp_bits)
	{
		return LEAF_NOT_SMOOTH;
	}
	if (!is_probable_prime128(p, 0) || !is_probable_prime128(q, 0))
	{
		return LEAF_NOT_SMOOTH;
	}
	lp[0] = (p < q) ? p : q;
	lp[1] = (p < q) ? q : p;
	return LEAF_SMOOTH;
}

//The cofactors that rho splits are listed (rho_list, rho_count) for k_rho, and are LEAF_TINYQS until then: rho takes long,
//and a few of them would hold up most blocks.
__global__ void k_classify(const uint64_t* __restrict__ c, size_t stride, const uint32_t* __restrict__ size, size_t count, Leaf_check check,
	unsigned char* __restrict__ result, uint64_t* __restrict__ lp, uint32_t* __restrict__ rho_list, uint32_t* __restrict__ rho_count)
{
	for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < count; i += (size_t)gridDim.x * blockDim.x)
	{
		const uint32_t n = size[i];
		const uint64_t* x = c + i * stride;
		const int bits = (n == 0) ? 1 : 64 * ((int)n - 1) + bits64(x[n - 1]);
		const uint64_t lo = (n > 0) ? x[0] : 0, hi = (n > 1) ? x[1] : 0;
		uint64_t primes[2] = { 0, 0 };
		unsigned char r;
		if (bits > check.max_cofactor_bits)
		{
			//Can still have three large primes, if composite and at least fblim^3.
			const bool at_least_min = (n > 2) || hi > check.three_lp_min[1] || (hi == check.three_lp_min[1] && lo >= check.three_lp_min[0]);
			if (bits <= check.three_lp_max_bits && at_least_min)
			{
				r = (n > 2) ? LEAF_CHECK_ON_CPU : (is_probable_prime128(lo, hi) ? LEAF_NOT_SMOOTH : LEAF_CANDIDATE);
			}
			else
			{
				r = LEAF_NOT_SMOOTH;
			}
		}
		else if (bits > check.two_lp_min_bits)
		{
			if (n > 2)
			{
				r = LEAF_CHECK_ON_CPU;
			}
			else if (is_probable_prime128(lo, hi))
			{
				r = LEAF_NOT_SMOOTH;
			}
			else
			{
				r = LEAF_TINYQS;
				if (bits <= 63 && (lo & 1))
				{
					rho_list[atomicAdd(rho_count, 1u)] = (uint32_t)i;
				}
			}
		}
		else if (bits > check.one_lp_bits)
		{
			r = LEAF_NOT_SMOOTH;
		}
		else if (n == 1 && lo == 1)
		{
			r = LEAF_SMOOTH;
		}
		else if (n <= 2 && (hi != 0 || lo > check.max_batch_prime) && is_probable_prime128(lo, hi))
		{
			primes[0] = lo; //bits <= one_lp_bits, so it fits in one limb
			r = LEAF_SMOOTH;
		}
		else
		{
			r = LEAF_CHECK_ON_CPU; //factored on the CPU, rarely needed
		}
		result[i] = r;
		lp[2 * i] = primes[0];
		lp[2 * i + 1] = primes[1];
	}
}

//Pollard rho on the cofactors that k_classify listed, one thread each, in small blocks so that they are spread over the GPU.
//A cofactor that rho does not split stays LEAF_TINYQS, as on the CPU.
__global__ void k_rho(const uint64_t* __restrict__ c, size_t stride, const uint32_t* __restrict__ rho_list, const uint32_t* __restrict__ rho_count, Leaf_check check,
	unsigned char* __restrict__ result, uint64_t* __restrict__ lp)
{
	const uint32_t count = *rho_count;
	for (uint32_t j = blockIdx.x * blockDim.x + threadIdx.x; j < count; j += gridDim.x * blockDim.x)
	{
		const uint32_t i = rho_list[j];
		const uint64_t n = c[(size_t)i * stride];
		const uint64_t factor = rho64(n);
		if (factor != 0)
		{
			uint64_t primes[2] = { 0, 0 };
			result[i] = check_two_primes(factor, n / factor, check, primes);
			lp[2 * (size_t)i] = primes[0];
			lp[2 * (size_t)i + 1] = primes[1];
		}
	}
}

//Factoring of the smooth leaves for the output, one warp per leaf: the value of the leaf without its large primes is factored as by
//factor_general of Smooth_factor.cpp: 2s and trial division (the small primes up to 128 bits, as factor_smooth, and below 2^16
//above), then rho until the parts are primes, and primes below 2^32 by is_prime32. But each rho runs 32 walks at once, one per
//lane with constants 1 to 32, and takes the factor of the first lane that finds one. The shortest of 32 walks is about sqrt(32)
//times shorter than one, and the lanes do the same steps. The factors found are primes below 2^32, so the factorization is the
//one that factor_general finds. A leaf that this does not factor (more than FACTOR_LIMBS limbs after its large primes, more than
//FACTOR_MAX factors, or walks of more than FACTOR_RHO_LIMIT) is factored on the CPU.
#define FACTOR_MAX 64 //the most factors of a value factored on the GPU
#define FACTOR_LIMBS 4 //the most limbs of a value factored on the GPU, and of the leaves
#define FACTOR_PENDING 24 //cofactors waiting to be split. After trial division, a value below 2^256 has at most 16 prime factors,
                          //and one below 2^128 at most 17, so this is never reached.
#define FACTOR_TRIAL_PRIMES 6541 //the odd primes below 2^16
#define FACTOR_POOL_PER_LEAF 8 //the room for factors in a batch, per leaf at first (usually a few percent of the leaves are smooth)
#define FACTOR_RHO_LIMIT (1 << 15) //the longest cycle length r of the walks before giving up, about 2^16 steps
#define FULL_WARP 0xFFFFFFFFu

//The factor of the lowest lane of the warp whose walk found one (found: the lanes that did), the same in all lanes.
__device__ __forceinline__ uint64_t first_found(unsigned int found, uint64_t factor)
{
	return __shfl_sync(FULL_WARP, factor, __ffs(found) - 1);
}

//A nontrivial factor of the odd composite n < 2^64, by Pollard-Brent rho with a walk per lane, or 0 if none was found with cycle
//lengths up to limit. Called by the whole warp, with the same n.
__device__ uint64_t warp_rho64(uint64_t n, uint64_t limit)
{
	const uint64_t ninv = montgomery_inverse64(n);
	uint64_t c = (threadIdx.x & 31) + 1;
	uint64_t y = 2, x = y, q = 1;
	for (uint64_t r = 1; r <= limit; r *= 2)
	{
		x = y;
		for (uint64_t i = 0; i < r; i++)
		{
			y = rho_step64(y, c, n, ninv);
		}
		for (uint64_t k = 0; k < r; k += RHO_GCD_INTERVAL)
		{
			const uint64_t saved_y = y;
			const uint64_t steps = (r - k < RHO_GCD_INTERVAL) ? r - k : RHO_GCD_INTERVAL;
			for (uint64_t i = 0; i < steps; i++)
			{
				y = rho_step64(y, c, n, ninv);
				q = montgomery_multiply64(q, (x > y) ? x - y : y - x, n, ninv);
			}
			uint64_t g = gcd64(q, n);
			if (g == n || g == 0)
			{
				//Several factors were found at once, go back and take one step at a time.
				uint64_t t = saved_y;
				g = 1;
				for (uint64_t i = 0; i < steps && g == 1; i++)
				{
					t = rho_step64(t, c, n, ninv);
					g = gcd64((x > t) ? x - t : t - x, n);
				}
			}
			const unsigned int found = __ballot_sync(FULL_WARP, g != 1 && g != n);
			if (found != 0)
			{
				return first_found(found, g);
			}
			if (g != 1)
			{
				//This walk found all factors at once: a new one, with another constant.
				c += 32;
				x = y = 2;
				q = 1;
			}
		}
	}
	return 0;
}

//x / d in place, for odd d that divides x (n limbs), by Hensel division. Returns the size of the quotient.
__device__ int divide_exact(uint64_t* x, int n, uint64_t d)
{
	const uint64_t inverse = 0 - montgomery_inverse64(d); //1 / d mod 2^64
	uint64_t borrow = 0;
	for (int k = 0; k < n; k++)
	{
		const uint64_t a = x[k];
		const uint64_t s = a - borrow;
		const uint64_t q = s * inverse;
		x[k] = q;
		borrow = __umul64hi(q, d) + (a < borrow); //q d = s + 2^64 umulhi(q, d), as q d = s mod 2^64
	}
	while (n > 0 && x[n - 1] == 0)
	{
		n--;
	}
	return n;
}

//Strong probable prime test to base b, for odd n < 2^32 that b < 2^32 does not divide, with Montgomery arithmetic.
__device__ bool sprp32(uint64_t n, uint64_t b)
{
	const uint64_t ninv = montgomery_inverse64(n);
	const uint64_t one = (0 - n) % n; //2^64 mod n
	const uint64_t minus_one = n - one;
	const uint64_t base = b * one % n; //less than 2^64, as both are below 2^32
	const int s = __ffsll((long long)(n - 1)) - 1;
	const uint64_t d = (n - 1) >> s;
	uint64_t x = one;
	for (int i = 63 - __clzll((long long)d); i >= 0; i--)
	{
		x = montgomery_multiply64(x, x, n, ninv);
		if ((d >> i) & 1)
		{
			x = montgomery_multiply64(x, base, n, ninv);
		}
	}
	if (x == one || x == minus_one)
	{
		return true;
	}
	for (int i = 1; i < s; i++)
	{
		x = montgomery_multiply64(x, x, n, ninv);
		if (x == minus_one)
		{
			return true;
		}
	}
	return false;
}

//is_prime32 (Miller-Rabin with the bases 2, 7 and 61, exact below 2^32), with Montgomery arithmetic, for odd n < 2^32.
__device__ bool is_odd_prime32(uint64_t n)
{
	if (n < 3)
	{
		return false;
	}
	const uint64_t bases[3] = { 2, 7, 61 };
	for (int b = 0; b < 3; b++)
	{
		if (n == bases[b])
		{
			return true;
		}
		if (n % bases[b] == 0 || !sprp32(n, bases[b]))
		{
			return false;
		}
	}
	return true;
}

//Numbers of up to FACTOR_LIMBS limbs, with their size n (the top limb nonzero, 0 limbs for 0), for the values of more than 128 bits.

//Whether odd p divides x, with inverse = 1 / p mod 2^64: the Hensel division of divide_exact leaves no borrow exactly then.
__device__ __forceinline__ bool divisible(const uint64_t* x, int n, uint64_t p, uint64_t inverse)
{
	uint64_t borrow = 0;
	for (int k = 0; k < n; k++)
	{
		const uint64_t a = x[k];
		const uint64_t q = (a - borrow) * inverse;
		borrow = __umul64hi(q, p) + (a < borrow);
	}
	return borrow == 0;
}

//c / d in place, for odd d of dn limbs that divides c (cn limbs), by Hensel division. Returns the size of the quotient.
__device__ int divide_exact_limbs(uint64_t* c, int cn, const uint64_t* d, int dn)
{
	if (dn == 1)
	{
		return divide_exact(c, cn, d[0]);
	}
	const uint64_t inverse = 0 - montgomery_inverse64(d[0]);
	const int qn = cn - dn + 1;
	uint64_t q[FACTOR_LIMBS];
	for (int i = 0; i < qn; i++)
	{
		const uint64_t qi = c[i] * inverse;
		q[i] = qi;
		//c -= qi d 2^(64 i)
		uint64_t borrow = 0;
		for (int j = 0; j < dn && i + j < cn; j++)
		{
			const uint64_t lo = qi * d[j], hi = __umul64hi(qi, d[j]);
			const uint64_t a = c[i + j];
			const uint64_t r = a - lo;
			const uint64_t b1 = (a < lo);
			c[i + j] = r - borrow;
			borrow = hi + b1 + (r < borrow);
		}
		for (int j = i + dn; j < cn; j++)
		{
			const uint64_t a = c[j];
			c[j] = a - borrow;
			borrow = (a < borrow);
		}
	}
	for (int i = 0; i < qn; i++)
	{
		c[i] = q[i];
	}
	int n = qn;
	while (n > 0 && c[n - 1] == 0)
	{
		n--;
	}
	return n;
}

//x >>= bits, for bits < 64 n, with the limbs above the new size zero. Returns the new size.
__device__ int shift_right_limbs(uint64_t* x, int n, int bits)
{
	const int words = bits / 64, b = bits % 64;
	for (int k = 0; k + words < n; k++)
	{
		const uint64_t low = x[k + words] >> b;
		const uint64_t high = (b != 0 && k + words + 1 < n) ? x[k + words + 1] << (64 - b) : 0;
		x[k] = low | high;
	}
	for (int k = (n - words > 0) ? n - words : 0; k < n; k++)
	{
		x[k] = 0;
	}
	n -= words;
	while (n > 0 && x[n - 1] == 0)
	{
		n--;
	}
	return n;
}

//The N-limb arithmetic of rho, modulo an odd n of N limbs (the top one nonzero), for numbers below n.

//hi:t = a b + t + carry, the new carry being hi.
__device__ __forceinline__ void multiply_add(uint64_t a, uint64_t b, uint64_t& t, uint64_t& carry)
{
	uint64_t lo = a * b, hi = __umul64hi(a, b);
	lo += t;
	hi += (lo < t);
	lo += carry;
	hi += (lo < carry);
	t = lo;
	carry = hi;
}

//x >= n
template <int N> __device__ __forceinline__ bool at_least(const uint64_t* x, const uint64_t* n)
{
	for (int k = N - 1; k > 0; k--)
	{
		if (x[k] != n[k])
		{
			return x[k] > n[k];
		}
	}
	return x[0] >= n[0];
}

//x -= n, returning the borrow.
template <int N> __device__ __forceinline__ uint64_t subtract_limbs(uint64_t* x, const uint64_t* n)
{
	uint64_t borrow = 0;
	for (int k = 0; k < N; k++)
	{
		const uint64_t a = x[k];
		const uint64_t r = a - n[k];
		const uint64_t b1 = (a < n[k]);
		x[k] = r - borrow;
		borrow = b1 + (r < borrow);
	}
	return borrow;
}

//a = a b / 2^(64 N) mod n, by Montgomery multiplication (CIOS).
template <int N> __device__ void montgomery_multiply_limbs(uint64_t* a, const uint64_t* b, const uint64_t* n, uint64_t ninv)
{
	uint64_t t[N + 1];
	for (int k = 0; k <= N; k++)
	{
		t[k] = 0;
	}
	for (int i = 0; i < N; i++)
	{
		//t += a b_i
		uint64_t carry = 0;
		for (int j = 0; j < N; j++)
		{
			multiply_add(a[j], b[i], t[j], carry);
		}
		t[N] += carry;
		uint64_t top = (t[N] < carry);
		//t = (t + m n) / 2^64, where m makes the lowest limb zero.
		const uint64_t m = t[0] * ninv;
		carry = 0;
		multiply_add(m, n[0], t[0], carry);
		for (int j = 1; j < N; j++)
		{
			multiply_add(m, n[j], t[j], carry);
			t[j - 1] = t[j];
		}
		t[N - 1] = t[N] + carry;
		t[N] = top + (t[N - 1] < carry);
	}
	if (t[N] != 0 || at_least<N>(t, n))
	{
		subtract_limbs<N>(t, n);
	}
	for (int k = 0; k < N; k++)
	{
		a[k] = t[k];
	}
}

//y = y^2 + c mod n
template <int N> __device__ __forceinline__ void rho_step_limbs(uint64_t* y, uint64_t c, const uint64_t* n, uint64_t ninv)
{
	montgomery_multiply_limbs<N>(y, y, n, ninv);
	uint64_t carry = c;
	for (int k = 0; k < N; k++)
	{
		y[k] += carry;
		carry = (y[k] < carry);
	}
	if (carry != 0 || at_least<N>(y, n))
	{
		subtract_limbs<N>(y, n);
	}
}

//d = |a - b|
template <int N> __device__ __forceinline__ void difference_limbs(const uint64_t* a, const uint64_t* b, uint64_t* d)
{
	for (int k = 0; k < N; k++)
	{
		d[k] = a[k];
	}
	if (subtract_limbs<N>(d, b) != 0)
	{
		for (int k = 0; k < N; k++)
		{
			d[k] = b[k];
		}
		subtract_limbs<N>(d, a);
	}
}

//b = gcd(a, b) for odd b, by the binary algorithm. a is changed.
template <int N> __device__ void gcd_limbs(uint64_t* a, uint64_t* b)
{
	for (;;)
	{
		int k = 0;
		while (k < N && a[k] == 0)
		{
			k++;
		}
		if (k == N)
		{
			return;
		}
		shift_right_limbs(a, N, 64 * k + __ffsll((long long)a[k]) - 1);
		if (!at_least<N>(a, b))
		{
			for (int i = 0; i < N; i++)
			{
				const uint64_t t = a[i];
				a[i] = b[i];
				b[i] = t;
			}
		}
		subtract_limbs<N>(a, b);
	}
}

template <int N> __device__ __forceinline__ bool equal_limbs(const uint64_t* a, const uint64_t* b)
{
	bool same = true;
	for (int k = 0; k < N; k++)
	{
		same &= (a[k] == b[k]);
	}
	return same;
}

template <int N> __device__ __forceinline__ bool is_one(const uint64_t* a)
{
	bool one = (a[0] == 1);
	for (int k = 1; k < N; k++)
	{
		one &= (a[k] == 0);
	}
	return one;
}

//warp_rho64 for an odd composite n of N limbs, with the factor in f (N limbs, zero above its size). Returns false if no factor
//was found.
template <int N> __device__ bool warp_rho_limbs(const uint64_t* n, uint64_t limit, uint64_t* f)
{
	const uint64_t ninv = montgomery_inverse64(n[0]);
	uint64_t c = (threadIdx.x & 31) + 1;
	uint64_t y[N], x[N], q[N], saved[N], d[N], g[N];
	for (int k = 0; k < N; k++)
	{
		y[k] = x[k] = q[k] = 0;
	}
	y[0] = x[0] = 2;
	q[0] = 1;
	for (uint64_t r = 1; r <= limit; r *= 2)
	{
		for (int k = 0; k < N; k++)
		{
			x[k] = y[k];
		}
		for (uint64_t i = 0; i < r; i++)
		{
			rho_step_limbs<N>(y, c, n, ninv);
		}
		for (uint64_t s = 0; s < r; s += RHO_GCD_INTERVAL)
		{
			for (int k = 0; k < N; k++)
			{
				saved[k] = y[k];
			}
			const uint64_t steps = (r - s < RHO_GCD_INTERVAL) ? r - s : RHO_GCD_INTERVAL;
			for (uint64_t i = 0; i < steps; i++)
			{
				rho_step_limbs<N>(y, c, n, ninv);
				difference_limbs<N>(x, y, d);
				montgomery_multiply_limbs<N>(q, d, n, ninv);
			}
			for (int k = 0; k < N; k++)
			{
				d[k] = q[k];
				g[k] = n[k];
			}
			gcd_limbs<N>(d, g);
			if (equal_limbs<N>(g, n))
			{
				//Several factors were found at once, go back and take one step at a time.
				for (uint64_t i = 0; i < steps && (i == 0 || is_one<N>(g)); i++)
				{
					rho_step_limbs<N>(saved, c, n, ninv);
					difference_limbs<N>(x, saved, d);
					for (int k = 0; k < N; k++)
					{
						g[k] = n[k];
					}
					gcd_limbs<N>(d, g);
				}
			}
			const bool one = is_one<N>(g), all = equal_limbs<N>(g, n);
			const unsigned int found = __ballot_sync(FULL_WARP, !one && !all);
			if (found != 0)
			{
				for (int k = 0; k < N; k++)
				{
					f[k] = first_found(found, g[k]);
				}
				return true;
			}
			if (!one)
			{
				//This walk found all factors at once: a new one, with another constant.
				c += 32;
				for (int k = 0; k < N; k++)
				{
					y[k] = x[k] = q[k] = 0;
				}
				y[0] = x[0] = 2;
				q[0] = 1;
			}
		}
	}
	return false;
}

//The odd primes below 2^16 for the trial division of the values of more than two limbs, and their inverses modulo 2^64.
__device__ uint64_t factor_trial[2 * FACTOR_TRIAL_PRIMES];

//The factors of the value x of n limbs (at most FACTOR_LIMBS, changed), by the whole warp (see above), in f by lane 0, in
//increasing order. Returns their number, or -1 if it was not factored.
__device__ int warp_factor(uint64_t* x, int n, uint32_t* f, uint64_t limit)
{
	const int lane = threadIdx.x & 31;
	if (n == 0)
	{
		return -1;
	}
	int count = 0;
	int twos = 0;
	while (x[twos / 64] == 0)
	{
		twos += 64;
	}
	twos += __ffsll((long long)x[twos / 64]) - 1;
	if (twos > FACTOR_MAX)
	{
		return -1;
	}
	for (int k = 0; k < twos; k++)
	{
		if (lane == 0)
		{
			f[count] = 2;
		}
		count++;
	}
	if (twos > 0)
	{
		n = shift_right_limbs(x, n, twos);
	}

	//Trial division by the first 32 odd primes (up to 137), a lane each, and above two limbs by the odd primes below 2^16 until the
	//rest has two limbs.
	for (int first = 0; first < FACTOR_TRIAL_PRIMES && (first == 0 || n > 2); first += 32)
	{
		const int k = first + lane;
		unsigned int divides = __ballot_sync(FULL_WARP, k < FACTOR_TRIAL_PRIMES && divisible(x, n, factor_trial[2 * k], factor_trial[2 * k + 1]));
		while (divides != 0)
		{
			const int i = first + __ffs(divides) - 1;
			divides &= divides - 1;
			const uint64_t p = factor_trial[2 * i], inverse = factor_trial[2 * i + 1];
			do
			{
				if (count >= FACTOR_MAX)
				{
					return -1;
				}
				if (lane == 0)
				{
					f[count] = (uint32_t)p;
				}
				count++;
				n = divide_exact(x, n, p);
			} while (divisible(x, n, p, inverse));
		}
	}

	//Split the rest. All cofactors are odd, and have no prime factors below 137.
	uint64_t pending[FACTOR_PENDING][FACTOR_LIMBS];
	int pending_size[FACTOR_PENDING];
	int num_pending = 0;
	if (!(n == 1 && x[0] == 1))
	{
		for (int k = 0; k < FACTOR_LIMBS; k++)
		{
			pending[0][k] = (k < n) ? x[k] : 0;
		}
		pending_size[0] = n;
		num_pending = 1;
	}
	while (num_pending > 0)
	{
		num_pending--;
		uint64_t c[FACTOR_LIMBS];
		for (int k = 0; k < FACTOR_LIMBS; k++)
		{
			c[k] = pending[num_pending][k];
		}
		const int cn = pending_size[num_pending];
		if (cn == 1 && c[0] < (1ULL << 32) && is_odd_prime32(c[0]))
		{
			if (count >= FACTOR_MAX)
			{
				return -1;
			}
			if (lane == 0)
			{
				f[count] = (uint32_t)c[0];
			}
			count++;
			continue;
		}
		if (num_pending + 2 > FACTOR_PENDING)
		{
			return -1;
		}
		uint64_t d[FACTOR_LIMBS] = { 0, 0, 0, 0 };
		bool found;
		if (cn == 1)
		{
			d[0] = warp_rho64(c[0], limit);
			found = (d[0] != 0);
		}
		else if (cn == 2)
		{
			found = warp_rho_limbs<2>(c, limit, d);
		}
		else if (cn == 3)
		{
			found = warp_rho_limbs<3>(c, limit, d);
		}
		else
		{
			found = warp_rho_limbs<4>(c, limit, d);
		}
		if (!found)
		{
			return -1;
		}
		int dn = cn;
		while (d[dn - 1] == 0)
		{
			dn--;
		}
		const int en = divide_exact_limbs(c, cn, d, dn);
		for (int k = 0; k < FACTOR_LIMBS; k++)
		{
			pending[num_pending][k] = d[k];
			pending[num_pending + 1][k] = (k < en) ? c[k] : 0;
		}
		pending_size[num_pending] = dn;
		pending_size[num_pending + 1] = en;
		num_pending += 2;
	}

	if (lane == 0)
	{
		for (int i = 1; i < count; i++)
		{
			const uint32_t y = f[i];
			int j = i - 1;
			while (j >= 0 && f[j] > y)
			{
				f[j + 1] = f[j];
				j--;
			}
			f[j + 1] = y;
		}
	}
	return count;
}

//The leaves to factor, in any order: LEAF_SMOOTH if they were checked (classes), and otherwise those with cofactor 1, except the
//unused leaves (1).
__global__ void k_list_smooth(const unsigned char* __restrict__ classes, const uint64_t* __restrict__ cofactors, size_t cofactor_stride,
	const uint32_t* __restrict__ cofactor_size, const uint64_t* __restrict__ leaves, size_t stride, const uint32_t* __restrict__ size, size_t count,
	uint32_t* __restrict__ list, uint32_t* __restrict__ list_count)
{
	for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < count; i += (size_t)gridDim.x * blockDim.x)
	{
		const bool smooth = (classes != NULL) ? classes[i] == LEAF_SMOOTH : (cofactor_size[i] == 1 && cofactors[i * cofactor_stride] == 1);
		if (smooth && !(size[i] == 1 && leaves[i * stride] == 1))
		{
			list[atomicAdd(list_count, 1u)] = (uint32_t)i;
		}
	}
}

#define FACTOR_BLOCK 128 //threads in a block of k_factor_smooth

//The factors of the listed leaves without their large primes (lp, two per leaf, or none), a warp each: factor_count[j] of them for list[j] (or -1 if it is
//left to the CPU), from factor_offset[j] in the pool, which has room for pool_capacity. pool_used counts the room they need, also
//beyond pool_capacity: those that do not fit are left to the CPU.
__global__ void __launch_bounds__(FACTOR_BLOCK) k_factor_smooth(const uint32_t* __restrict__ list, const uint32_t* __restrict__ list_count, const uint64_t* __restrict__ leaves,
	size_t stride, const uint32_t* __restrict__ size, const uint64_t* __restrict__ lp, int* __restrict__ factor_count, uint32_t* __restrict__ factor_offset,
	uint32_t* __restrict__ pool, uint32_t* __restrict__ pool_used, uint32_t pool_capacity)
{
	__shared__ uint32_t found[FACTOR_BLOCK / 32][FACTOR_MAX];
	uint32_t* f = found[threadIdx.x / 32];
	const uint32_t count = *list_count;
	const uint32_t warps = gridDim.x * (FACTOR_BLOCK / 32);
	for (uint32_t j = blockIdx.x * (FACTOR_BLOCK / 32) + threadIdx.x / 32; j < count; j += warps)
	{
		const uint32_t i = list[j];
		int n = (int)size[i];
		int result = -1;
		if (n <= FACTOR_LIMBS)
		{
			uint64_t x[FACTOR_LIMBS];
			for (int k = 0; k < FACTOR_LIMBS; k++)
			{
				x[k] = (k < n) ? leaves[(size_t)i * stride + k] : 0;
			}
			for (int k = 0; k < 2; k++)
			{
				const uint64_t p = (lp != NULL) ? lp[2 * (size_t)i + k] : 0;
				if (p > 1)
				{
					n = divide_exact(x, n, p);
				}
			}
			result = warp_factor(x, n, f, FACTOR_RHO_LIMIT);
		}
		if ((threadIdx.x & 31) == 0)
		{
			uint32_t offset = 0;
			if (result > 0)
			{
				offset = atomicAdd(pool_used, (uint32_t)result);
				if (offset + result > pool_capacity)
				{
					result = -1;
				}
				for (int k = 0; k < result; k++)
				{
					pool[offset + k] = f[k];
				}
			}
			factor_count[j] = result;
			factor_offset[j] = offset;
		}
		__syncwarp();
	}
}

#define OTHER_MAX_DEGREE 8 //the largest degree of the other side's polynomial on the GPU
#define OTHER_WIDTH 8 //limbs of the two's complement numbers in which the other side's norms are computed

//x *= u modulo 2^(64 OTHER_WIDTH), which is the two's complement product as long as it fits.
__device__ __forceinline__ void multiply_scalar_wide(uint64_t* x, uint64_t u)
{
	uint64_t carry = 0;
	for (int k = 0; k < OTHER_WIDTH; k++)
	{
		const uint64_t lo = x[k] * u, hi = __umul64hi(x[k], u);
		x[k] = lo + carry;
		carry = hi + (x[k] < lo);
	}
}

__device__ __forceinline__ void negate_wide(uint64_t* x)
{
	uint64_t carry = 1;
	for (int k = 0; k < OTHER_WIDTH; k++)
	{
		x[k] = ~x[k] + carry;
		carry = (x[k] < carry);
	}
}

//The norms |F(a, b)| = |sum c_k a^k b^(d - k)| of the other side of the listed leaves (ab: a and b of each leaf; coefficients: of
//OTHER_WIDTH limbs each, two's complement, from c_0), as numbers of up to FACTOR_LIMBS limbs (values, at the stride FACTOR_LIMBS, and
//sizes, by leaf). Norms that may not fit in OTHER_WIDTH limbs (from the sizes of a, b and the coefficients), or do not fit in
//FACTOR_LIMBS, get the size FACTOR_LIMBS + 1, which k_factor_smooth leaves to the CPU.
__global__ void k_other_norms(const uint32_t* __restrict__ list, const uint32_t* __restrict__ list_count, const int64_t* __restrict__ ab,
	const uint64_t* __restrict__ coefficients, int degree, int coefficient_bits, uint64_t* __restrict__ values, uint32_t* __restrict__ sizes)
{
	const uint32_t count = *list_count;
	for (uint32_t j = blockIdx.x * blockDim.x + threadIdx.x; j < count; j += gridDim.x * blockDim.x)
	{
		const uint32_t i = list[j];
		const int64_t a = ab[2 * (size_t)i], b = ab[2 * (size_t)i + 1];
		const uint64_t abs_a = (a < 0) ? 0 - (uint64_t)a : (uint64_t)a, abs_b = (b < 0) ? 0 - (uint64_t)b : (uint64_t)b;
		const uint64_t larger = (abs_a > abs_b) ? abs_a : abs_b;
		const int larger_bits = 64 - __clzll((long long)larger);
		if (coefficient_bits + degree * larger_bits + 4 > 64 * OTHER_WIDTH - 1)
		{
			sizes[i] = FACTOR_LIMBS + 1;
			continue;
		}
		uint64_t sum[OTHER_WIDTH], term[OTHER_WIDTH];
		for (int k = 0; k < OTHER_WIDTH; k++)
		{
			sum[k] = 0;
		}
		for (int e = 0; e <= degree; e++)
		{
			for (int k = 0; k < OTHER_WIDTH; k++)
			{
				term[k] = coefficients[e * OTHER_WIDTH + k];
			}
			for (int k = 0; k < e; k++)
			{
				multiply_scalar_wide(term, abs_a);
			}
			for (int k = e; k < degree; k++)
			{
				multiply_scalar_wide(term, abs_b);
			}
			if (((a < 0) && (e & 1)) != ((b < 0) && ((degree - e) & 1)))
			{
				negate_wide(term);
			}
			uint64_t carry = 0;
			for (int k = 0; k < OTHER_WIDTH; k++)
			{
				const uint64_t s = sum[k] + term[k];
				const uint64_t c1 = (s < term[k]);
				sum[k] = s + carry;
				carry = c1 + (sum[k] < s);
			}
		}
		if (sum[OTHER_WIDTH - 1] >> 63)
		{
			negate_wide(sum);
		}
		int n = OTHER_WIDTH;
		while (n > 0 && sum[n - 1] == 0)
		{
			n--;
		}
		if (n > FACTOR_LIMBS)
		{
			sizes[i] = FACTOR_LIMBS + 1;
			continue;
		}
		for (int k = 0; k < FACTOR_LIMBS; k++)
		{
			values[(size_t)i * FACTOR_LIMBS + k] = sum[k];
		}
		sizes[i] = (uint32_t)n;
	}
}

//The cofactors that the CPU still needs after the checks (LEAF_CANDIDATE, LEAF_TINYQS and LEAF_CHECK_ON_CPU), gathered for the
//copy to the host, in any order: their leaves (list, after their count), sizes and limbs (K each).
__global__ void k_gather_cofactors(const unsigned char* __restrict__ classes, const uint64_t* __restrict__ cofactors, const uint32_t* __restrict__ size, int K,
	size_t count, uint32_t* __restrict__ list, uint32_t* __restrict__ list_count, uint32_t* __restrict__ gathered_size, uint64_t* __restrict__ gathered)
{
	for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < count; i += (size_t)gridDim.x * blockDim.x)
	{
		const unsigned char c = classes[i];
		if (c == LEAF_CANDIDATE || c == LEAF_TINYQS || c == LEAF_CHECK_ON_CPU)
		{
			const uint32_t j = atomicAdd(list_count, 1u);
			list[j] = (uint32_t)i;
			gathered_size[j] = size[i];
			for (int k = 0; k < K; k++)
			{
				gathered[(size_t)j * K + k] = cofactors[i * K + k];
			}
		}
	}
}

} //namespace

class Gpu_tree_backend : public Cpu_tree_backend
{
public:
	Gpu_tree_backend(int levels, int max_bits_in_input);
	~Gpu_tree_backend();
	mp_bitcnt_t Build_product_tree();
	void Scaled_root(mpz_srcptr prime_product);
	void Remainder_tree();
	void Leaf_cofactors();
	bool Set_leaf_check(const Leaf_check& check);
	const unsigned char* Leaf_classes();
	const unsigned long long* Leaf_large_primes();
	bool Set_leaf_factoring();
	bool Leaf_factors_ready();
	int Leaf_factors(int index, const unsigned int** factors);

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
	bool checking; //the cofactors, with check
	bool factoring; //the smooth leaves, with Start_factoring
	Leaf_check check;
	bool classes_valid; //for this batch
	unsigned char* classes; //on the GPU, and host_classes
	size_t classes_capacity;
	uint64_t* large_primes; //on the GPU, two per leaf, and host_large_primes
	size_t large_primes_capacity;
	uint32_t* rho_list; //the leaves for k_rho, and their count at rho_list[0]
	size_t rho_list_capacity;
	cudaStream_t rho_stream; //k_rho runs on a stream of its own, after classified (on the stream of the thread)
	cudaEvent_t classified, rho_done; //k_rho of the batch, which the thread that checks the batch waits for (Leaf_classes)
	bool rho_pending;
	void Finish_rho();
	//The factoring of the smooth leaves, after k_rho on rho_stream: the leaves (factor_list, after their count and the room their
	//factors need in the pool), and for each its number of factors and their offset in the pool (factor_primes).
	uint32_t* factor_list;
	int* factor_counts;
	uint32_t *factor_offsets, *factor_primes;
	size_t factor_list_capacity, factor_counts_capacity, factor_offsets_capacity, factor_primes_capacity;
	uint32_t factor_pool; //the room in the pool for this batch: FACTOR_POOL_PER_LEAF per leaf, or what the last batch needed
	cudaEvent_t factored;
	bool factors_pending, factors_valid;
	void Start_factoring(int K);
	std::vector<uint32_t> host_factor_list, host_factor_offsets, host_factor_primes;
	std::vector<int> host_factor_counts, host_factor_entry; //host_factor_entry: the entry of each leaf in the list, or -1
	Pinned<unsigned char> host_classes;
	Pinned<unsigned long long> host_large_primes;
	//The cofactors that the CPU needs, gathered on the GPU (k_gather_cofactors): their leaves (after their count), sizes and limbs.
	uint32_t *gather_list, *gather_sizes;
	uint64_t* gather_limbs;
	size_t gather_list_capacity, gather_sizes_capacity, gather_limbs_capacity;
	Pinned<uint32_t> host_gather_list;
	std::vector<Depth> depth_data;
	uint32_t* parent_transform;
	size_t parent_transform_capacity;
	uint32_t* pair_transform;
	size_t pair_transform_capacity;
	Gpu_level scratch_level; //cyclic products, and exact truncated products
	size_t scratch_limb_capacity, scratch_size_capacity;
	uint32_t *shift, *precision;
	size_t shift_capacity, precision_capacity;
	Pinned<uint64_t> host_limbs;
	Pinned<uint32_t> host_sizes;
	std::vector<uint32_t> host_words;
};

//The table of factor_trial, once.
static void set_trial_primes()
{
	static std::once_flag done;
	std::call_once(done, []
	{
		const uint32_t limit = 1 << 16;
		std::vector<bool> composite(limit, false);
		std::vector<uint64_t> table;
		for (uint32_t i = 3; i < limit; i += 2)
		{
			if (composite[i])
			{
				continue;
			}
			uint64_t inverse = i; //1 / i mod 2^64, by Newton iteration from 3 bits
			for (int k = 0; k < 5; k++)
			{
				inverse *= 2 - i * inverse;
			}
			table.push_back(i);
			table.push_back(inverse);
			for (uint32_t j = i * i; j < limit; j += 2 * i)
			{
				composite[j] = true;
			}
		}
		if (table.size() != 2 * FACTOR_TRIAL_PRIMES)
		{
			fprintf(stderr, "Gpu_tree_backend: %zu odd primes below 2^16\n", table.size() / 2);
			exit(3);
		}
		CUDA_CHECK(cudaMemcpyToSymbol(factor_trial, table.data(), table.size() * sizeof(uint64_t)));
	});
}

Gpu_tree_backend::Gpu_tree_backend(int levels, int max_bits_in_input) : Cpu_tree_backend(new Gpu_tree_products, levels, max_bits_in_input, false)
{
	gpu_start();
	set_trial_primes();
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
	checking = false;
	factoring = false;
	gather_list = gather_sizes = NULL;
	gather_limbs = NULL;
	gather_list_capacity = gather_sizes_capacity = gather_limbs_capacity = 0;
	classes_valid = false;
	classes = NULL;
	large_primes = NULL;
	rho_list = NULL;
	classes_capacity = large_primes_capacity = rho_list_capacity = 0;
	CUDA_CHECK(cudaStreamCreateWithFlags(&rho_stream, cudaStreamNonBlocking));
	CUDA_CHECK(cudaEventCreateWithFlags(&classified, cudaEventDisableTiming));
	CUDA_CHECK(cudaEventCreateWithFlags(&rho_done, cudaEventDisableTiming));
	rho_pending = false;
	factor_list = NULL;
	factor_counts = NULL;
	factor_offsets = factor_primes = NULL;
	factor_list_capacity = factor_counts_capacity = factor_offsets_capacity = factor_primes_capacity = 0;
	factor_pool = 0;
	CUDA_CHECK(cudaEventCreateWithFlags(&factored, cudaEventDisableTiming));
	factors_pending = factors_valid = false;
}

Gpu_tree_backend::~Gpu_tree_backend()
{
	//The kernels in the background read the leaves and the cofactors: wait for them before anything is freed.
	if (rho_pending)
	{
		cudaEventSynchronize(rho_done);
	}
	if (factors_pending)
	{
		cudaEventSynchronize(factored);
	}
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
	cudaEventDestroy(factored);
	cudaFree(gather_list);
	cudaFree(gather_sizes);
	cudaFree(gather_limbs);
	cudaFree(factor_list);
	cudaFree(factor_counts);
	cudaFree(factor_offsets);
	cudaFree(factor_primes);
	cudaEventDestroy(rho_done);
	cudaEventDestroy(classified);
	cudaStreamDestroy(rho_stream);
	cudaFree(classes);
	cudaFree(large_primes);
	cudaFree(rho_list);
}

bool Gpu_tree_backend::Set_leaf_check(const Leaf_check& leaf_check)
{
	check = leaf_check;
	checking = true;
	return true;
}

//The results of k_rho, in the thread that checks the batch: by then it is usually done.
void Gpu_tree_backend::Finish_rho()
{
	if (!rho_pending)
	{
		return;
	}
	gpu_start();
	CUDA_CHECK(cudaEventSynchronize(rho_done));
	host_classes.resize(num_leaves);
	host_large_primes.resize(2 * (size_t)num_leaves);
	CUDA_CHECK(cudaMemcpy(host_classes.data(), classes, num_leaves, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaMemcpy(host_large_primes.data(), large_primes, 2 * (size_t)num_leaves * sizeof(uint64_t), cudaMemcpyDeviceToHost));
	rho_pending = false;
}

const unsigned char* Gpu_tree_backend::Leaf_classes()
{
	Finish_rho();
	return classes_valid ? host_classes.data() : NULL;
}

const unsigned long long* Gpu_tree_backend::Leaf_large_primes()
{
	Finish_rho();
	return classes_valid ? host_large_primes.data() : NULL;
}

bool Gpu_tree_backend::Set_leaf_factoring()
{
	factoring = true;
	return true;
}

//The factoring of the smooth leaves on rho_stream, after what is already on it (k_rho) and on the stream of the thread (the
//cofactors), a warp per leaf.
void Gpu_tree_backend::Start_factoring(int K)
{
	Depth& leaves = depth_data[levels - 1];
	CUDA_CHECK(cudaEventRecord(classified, 0));
	CUDA_CHECK(cudaStreamWaitEvent(rho_stream, classified, 0));
	factor_pool = (factor_pool > FACTOR_POOL_PER_LEAF * (uint32_t)num_leaves) ? factor_pool : FACTOR_POOL_PER_LEAF * (uint32_t)num_leaves;
	reserve(factor_list, factor_list_capacity, (size_t)num_leaves + 2);
	reserve(factor_counts, factor_counts_capacity, (size_t)num_leaves);
	reserve(factor_offsets, factor_offsets_capacity, (size_t)num_leaves);
	reserve(factor_primes, factor_primes_capacity, (size_t)factor_pool);
	CUDA_CHECK(cudaMemsetAsync(factor_list, 0, 2 * sizeof(uint32_t), rho_stream));
	k_list_smooth<<<grid_size_leaves(num_leaves), 128, 0, rho_stream>>>(checking ? classes : NULL, cofactors.limbs, K, cofactors.size, leaves.level.limbs,
		leaves.level.stride, leaves.level.size, num_leaves, factor_list + 2, factor_list);
	k_factor_smooth<<<1024, FACTOR_BLOCK, 0, rho_stream>>>(factor_list + 2, factor_list, leaves.level.limbs, leaves.level.stride, leaves.level.size,
		checking ? large_primes : NULL, factor_counts, factor_offsets, factor_primes, factor_list + 1, factor_pool);
	CUDA_CHECK(cudaGetLastError());
	CUDA_CHECK(cudaEventRecord(factored, rho_stream));
	factors_pending = true;
}

//The factors of the smooth leaves, in the thread that checks the batch, which waits for them when it writes the relations.
bool Gpu_tree_backend::Leaf_factors_ready()
{
	if (!factors_pending)
	{
		return factors_valid;
	}
	gpu_start();
	CUDA_CHECK(cudaEventSynchronize(factored));
	uint32_t header[2]; //the number of leaves listed, and the room their factors need
	CUDA_CHECK(cudaMemcpy(header, factor_list, sizeof(header), cudaMemcpyDeviceToHost));
	const uint32_t count = header[0];
	host_factor_list.resize(count);
	host_factor_counts.resize(count);
	host_factor_offsets.resize(count);
	host_factor_primes.resize((header[1] < factor_pool) ? header[1] : factor_pool);
	if (count > 0)
	{
		CUDA_CHECK(cudaMemcpy(host_factor_list.data(), factor_list + 2, count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(host_factor_counts.data(), factor_counts, count * sizeof(int), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(host_factor_offsets.data(), factor_offsets, count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(host_factor_primes.data(), factor_primes, host_factor_primes.size() * sizeof(uint32_t), cudaMemcpyDeviceToHost));
	}
	if (header[1] > factor_pool)
	{
		//The pool was too small, and the CPU factors the leaves that did not fit. Room for all of them in the next batch.
		factor_pool = header[1] + header[1] / 4;
	}
	host_factor_entry.assign(num_leaves, -1);
	for (uint32_t j = 0; j < count; j++)
	{
		host_factor_entry[host_factor_list[j]] = (int)j;
	}
	factors_pending = false;
	factors_valid = true;
	return true;
}

int Gpu_tree_backend::Leaf_factors(int index, const unsigned int** factors)
{
	const int j = factors_valid ? host_factor_entry[index] : -1;
	if (j < 0)
	{
		return -1;
	}
	*factors = host_factor_primes.data() + host_factor_offsets[j];
	return host_factor_counts[j];
}

void Gpu_tree_backend::Set_stride(Depth& d, size_t count, size_t stride)
{
	d.level.stride = stride;
	d.level.count = count;
	reserve(d.level.limbs, d.limb_capacity, count * stride);
}

//Copy the nodes of a depth from the CPU to the GPU.
//The numbers are copied as wide as the largest of them (the stride can have room to spare, as for the leaves), into the scratch
//level, and spread out to the stride on the GPU.
void Gpu_tree_backend::Upload(int depth)
{
	Depth& d = depth_data[depth];
	const int first = (1 << depth) - 1;
	const size_t count = d.level.count, stride = d.level.stride;
	size_t width = 1;
	for (size_t r = 0; r < count; r++)
	{
		const size_t n = mpz_size(node[first + r].get_mpz_t());
		width = (n > width) ? n : width;
	}
	if (width > stride)
	{
		fprintf(stderr, "Gpu_tree_backend: a node of %zu limbs does not fit in %zu\n", width, stride);
		exit(3);
	}
	host_limbs.resize(count * width);
	host_sizes.resize(count);
	for (size_t r = 0; r < count; r++)
	{
		mpz_srcptr a = node[first + r].get_mpz_t();
		const size_t n = mpz_size(a);
		uint64_t* x = host_limbs.data() + r * width;
		memcpy(x, a->_mp_d, n * sizeof(uint64_t));
		memset(x + n, 0, (width - n) * sizeof(uint64_t));
		host_sizes[r] = (uint32_t)n;
	}
	reserve(scratch_level.limbs, scratch_limb_capacity, count * width);
	CUDA_CHECK(cudaMemsetAsync(d.level.limbs, 0, count * stride * sizeof(uint64_t)));
	CUDA_CHECK(cudaMemcpy(scratch_level.limbs, host_limbs.data(), host_limbs.size() * sizeof(uint64_t), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy2DAsync(d.level.limbs, stride * sizeof(uint64_t), scratch_level.limbs, width * sizeof(uint64_t), width * sizeof(uint64_t), count,
		cudaMemcpyDeviceToDevice));
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
	gpu_start();
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
			gpu->Combine_products(p, PRODUCTS_ADJACENT, NULL, X, parent_transform, count, false, parents.level, 0);
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
	gpu_start();
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

__global__ void k_shifts(const uint32_t* __restrict__ parent_bits, const uint32_t* __restrict__ node_bits, size_t count, uint32_t guard, uint32_t* __restrict__ shift,
	uint32_t* __restrict__ precision)
{
	for (size_t r = (size_t)blockIdx.x * blockDim.x + threadIdx.x; r < count; r += (size_t)gridDim.x * blockDim.x)
	{
		shift[r] = parent_bits[r / 2] - node_bits[r];
		precision[r] = parent_bits[r / 2] + guard;
	}
}

//The shifts and precisions of the middle products of the nodes of a depth: precision = parent bits + guard, and
//shift = parent bits - node bits, from the sizes of the product tree on the GPU.
void Gpu_tree_backend::Shifts(int depth)
{
	const size_t count = (size_t)1 << depth;
	reserve(shift, shift_capacity, count);
	reserve(precision, precision_capacity, count);
	const unsigned int blocks = (unsigned int)((count + 255) / 256 < 65536 ? (count + 255) / 256 : 65536);
	k_shifts<<<blocks, 256>>>(depth_data[depth - 1].bits, depth_data[depth].bits, count, (uint32_t)guard, shift, precision);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_tree_backend::Remainder_tree()
{
	gpu_start();
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
			gpu->Combine_products(p, PRODUCTS_SWAP, parent_transform, nodes.transform, pair_transform, count, true, scratch_level, 0);
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
	gpu_start();
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
			classes_valid = checking;
			factors_pending = factors_valid = false;
			if (checking)
			{
				reserve(classes, classes_capacity, (size_t)num_leaves);
				reserve(large_primes, large_primes_capacity, 2 * (size_t)num_leaves);
				reserve(rho_list, rho_list_capacity, (size_t)num_leaves + 1);
				CUDA_CHECK(cudaMemsetAsync(rho_list, 0, sizeof(uint32_t)));
				k_classify<<<grid_size_leaves(num_leaves), 128>>>(cofactors.limbs, K, cofactors.size, num_leaves, check, classes, large_primes, rho_list + 1, rho_list);
				CUDA_CHECK(cudaGetLastError());
			}
			//The cofactors as numbers: only those that the CPU still needs, if they were checked (the classes and large primes are
			//copied after k_rho, by Finish_rho); otherwise all.
			size_t count = num_leaves;
			const uint32_t* leaf = NULL;
			if (checking)
			{
				reserve(gather_list, gather_list_capacity, (size_t)num_leaves + 1);
				reserve(gather_sizes, gather_sizes_capacity, (size_t)num_leaves);
				reserve(gather_limbs, gather_limbs_capacity, (size_t)num_leaves * K);
				CUDA_CHECK(cudaMemsetAsync(gather_list, 0, sizeof(uint32_t)));
				k_gather_cofactors<<<grid_size_leaves(num_leaves), 128>>>(classes, cofactors.limbs, cofactors.size, K, num_leaves, gather_list + 1, gather_list,
					gather_sizes, gather_limbs);
				CUDA_CHECK(cudaGetLastError());
				uint32_t gathered = 0;
				CUDA_CHECK(cudaMemcpy(&gathered, gather_list, sizeof(uint32_t), cudaMemcpyDeviceToHost));
				count = gathered;
				host_gather_list.resize(count);
				if (count > 0)
				{
					CUDA_CHECK(cudaMemcpy(host_gather_list.data(), gather_list + 1, count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
				}
				leaf = host_gather_list.data();
			}
			host_limbs.resize(count * K);
			host_sizes.resize(count);
			if (count > 0)
			{
				CUDA_CHECK(cudaMemcpy(host_limbs.data(), checking ? gather_limbs : cofactors.limbs, count * K * sizeof(uint64_t), cudaMemcpyDeviceToHost));
				CUDA_CHECK(cudaMemcpy(host_sizes.data(), checking ? gather_sizes : cofactors.size, count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
			}
			parallel_for((int)count, 4096, [&](int j)
			{
				mpz_ptr a = cofactor[(leaf != NULL) ? leaf[j] : j].get_mpz_t();
				const int n = (int)host_sizes[j];
				if (a->_mp_alloc < n)
				{
					mpz_realloc2(a, n * GMP_LIMB_BITS);
				}
				memcpy(a->_mp_d, host_limbs.data() + (size_t)j * K, n * sizeof(uint64_t));
				a->_mp_size = n;
			});
			if (checking)
			{
				//Rho on the cofactors that need it, which is slow on the GPU (long chains of dependent products) but uses few of its
				//multiprocessors: the thread that checks the batch waits for it, so that this one can go on with the next batch. Until
				//then those cofactors are LEAF_TINYQS, and they are among the ones copied above.
				CUDA_CHECK(cudaEventRecord(classified, 0));
				CUDA_CHECK(cudaStreamWaitEvent(rho_stream, classified, 0));
				k_rho<<<256, 32, 0, rho_stream>>>(cofactors.limbs, K, rho_list + 1, rho_list, check, classes, large_primes);
				CUDA_CHECK(cudaGetLastError());
				CUDA_CHECK(cudaEventRecord(rho_done, rho_stream));
				rho_pending = true;
			}
			if (factoring)
			{
				//Then the factoring of the smooth leaves for the output, after rho on the same stream.
				Start_factoring(K);
			}
			return;
		}
	}
	//Large leaves: on the CPU, from the fractions of their parents. Their cofactors are not checked here.
	classes_valid = false;
	factors_pending = factors_valid = false;
	Download(levels - 2);
	Cpu_tree_backend::Leaf_cofactors();
}

#define OTHER_JOB_CAPACITY 4096 //norms of a job for which the memory is allocated up front (a batch of 2^19 relations has about 1500)

//The factoring of the other side's norms of the relations found, as jobs apart from the trees: a job (the relations found in a batch)
//is copied to the GPU, and its norms computed and factored (k_other_norms, k_factor_smooth, a warp per norm), on a stream of its own,
//so that the jobs of several batches run at once. A norm takes several splits by Pollard rho, which take
//tens of ms while the trees keep the GPU busy, so the results are only waited for when the relations are written, and the batches
//go on meanwhile. Submit is called by one thread, and Results by one thread, in the order of the jobs.
class Gpu_other_side_factoring : public Other_side_factoring
{
public:
	Gpu_other_side_factoring(const std::vector<uint64_t>& words, int degree, int coefficient_bits);
	~Gpu_other_side_factoring();
	int Submit(const long long* a, const long long* b, int count);
	int Submit_values(const unsigned long long* limbs, const unsigned int* sizes, int count);
	void Results(int job, std::vector<int>& counts, std::vector<unsigned int>& offsets, std::vector<unsigned int>& primes);

private:
	struct Job
	{
		cudaStream_t stream;
		cudaEvent_t done;
		bool busy; //submitted, and its results not taken yet
		int count;
		Pinned<int64_t> host_ab; //a and b of the norms
		Pinned<uint64_t> host_values; //or the values themselves (Submit_values)
		Pinned<uint32_t> host_sizes;
		Pinned<uint32_t> host_header; //the number of norms, and the room their factors need in the pool
		Pinned<int> host_counts;
		Pinned<uint32_t> host_offsets, host_primes;
		int64_t* ab;
		uint64_t* values;
		uint32_t *header, *list, *sizes, *offsets, *primes;
		int* counts;
		size_t capacity; //norms
	};
	void Reserve(Job& j, size_t count);
	int Launch(const long long* a, const long long* b, const unsigned long long* limbs, const unsigned int* sizes, int count);
	Job jobs[OTHER_JOBS];
	long long submitted, taken;
	std::mutex lock;
	std::condition_variable changed;
	uint64_t* coefficients;
	int degree, coefficient_bits;
};

Gpu_other_side_factoring::Gpu_other_side_factoring(const std::vector<uint64_t>& words, int degree, int coefficient_bits)
	: submitted(0), taken(0), degree(degree), coefficient_bits(coefficient_bits)
{
	gpu_start();
	set_trial_primes();
	CUDA_CHECK(cudaMalloc(&coefficients, words.size() * sizeof(uint64_t)));
	CUDA_CHECK(cudaMemcpy(coefficients, words.data(), words.size() * sizeof(uint64_t), cudaMemcpyHostToDevice));
	//At the highest priority, which was a few percent faster: the jobs are small, and their results are waited for.
	int least_priority = 0, greatest_priority = 0;
	CUDA_CHECK(cudaDeviceGetStreamPriorityRange(&least_priority, &greatest_priority));
	for (Job& j : jobs)
	{
		CUDA_CHECK(cudaStreamCreateWithPriority(&j.stream, cudaStreamNonBlocking, greatest_priority));
		CUDA_CHECK(cudaEventCreateWithFlags(&j.done, cudaEventDisableTiming));
		j.busy = false;
		j.count = 0;
		j.ab = NULL;
		j.values = NULL;
		j.header = j.list = j.sizes = j.offsets = j.primes = NULL;
		j.counts = NULL;
		j.capacity = 0;
		CUDA_CHECK(cudaMalloc(&j.header, 2 * sizeof(uint32_t)));
		j.host_header.resize(2);
		Reserve(j, OTHER_JOB_CAPACITY);
	}
	//The coefficients too are on the GPU before the first job.
	CUDA_CHECK(cudaDeviceSynchronize());
}

//Room for count norms in a job, their factors, and the list of them (0, 1, ...), which k_other_norms and k_factor_smooth take.
//The memory is allocated up front, so that it is not freed and allocated while the GPU is busy.
void Gpu_other_side_factoring::Reserve(Job& j, size_t count)
{
	if (count <= j.capacity)
	{
		return;
	}
	cudaFree(j.ab);
	cudaFree(j.values);
	cudaFree(j.list);
	cudaFree(j.sizes);
	cudaFree(j.offsets);
	cudaFree(j.primes);
	cudaFree(j.counts);
	CUDA_CHECK(cudaMalloc(&j.ab, 2 * count * sizeof(int64_t)));
	CUDA_CHECK(cudaMalloc(&j.values, count * FACTOR_LIMBS * sizeof(uint64_t)));
	CUDA_CHECK(cudaMalloc(&j.list, count * sizeof(uint32_t)));
	CUDA_CHECK(cudaMalloc(&j.sizes, count * sizeof(uint32_t)));
	CUDA_CHECK(cudaMalloc(&j.offsets, count * sizeof(uint32_t)));
	CUDA_CHECK(cudaMalloc(&j.primes, count * FACTOR_MAX * sizeof(uint32_t)));
	CUDA_CHECK(cudaMalloc(&j.counts, count * sizeof(int)));
	std::vector<uint32_t> identity(count);
	for (size_t i = 0; i < count; i++)
	{
		identity[i] = (uint32_t)i;
	}
	//On the stream of the job, and waited for: cudaMemcpy from memory that is not page-locked can return before the data is on the
	//GPU, and the kernels of the job (on its own stream) then read a list that is not there yet.
	CUDA_CHECK(cudaMemcpyAsync(j.list, identity.data(), count * sizeof(uint32_t), cudaMemcpyHostToDevice, j.stream));
	CUDA_CHECK(cudaStreamSynchronize(j.stream));
	j.host_ab.resize(2 * count);
	j.host_values.resize(count * FACTOR_LIMBS);
	j.host_sizes.resize(count);
	j.host_counts.resize(count);
	j.host_offsets.resize(count);
	j.host_primes.resize(count * FACTOR_MAX);
	j.capacity = count;
}

Gpu_other_side_factoring::~Gpu_other_side_factoring()
{
	for (Job& j : jobs)
	{
		cudaStreamSynchronize(j.stream);
		cudaFree(j.ab);
		cudaFree(j.values);
		cudaFree(j.header);
		cudaFree(j.list);
		cudaFree(j.sizes);
		cudaFree(j.offsets);
		cudaFree(j.primes);
		cudaFree(j.counts);
		cudaEventDestroy(j.done);
		cudaStreamDestroy(j.stream);
	}
	cudaFree(coefficients);
}

//Start a job: the norms of the relations (a, b) and their factors, at most FACTOR_MAX each. Waits for the oldest job's results to be
//taken if all OTHER_JOBS are in use. Returns the number of the job.
int Gpu_other_side_factoring::Submit(const long long* a, const long long* b, int count)
{
	return Launch(a, b, NULL, NULL, count);
}

int Gpu_other_side_factoring::Submit_values(const unsigned long long* limbs, const unsigned int* sizes, int count)
{
	return Launch(NULL, NULL, limbs, sizes, count);
}

//A job of the norms of (a, b), or (if a is NULL) of the values given.
int Gpu_other_side_factoring::Launch(const long long* a, const long long* b, const unsigned long long* limbs, const unsigned int* sizes, int count)
{
	std::unique_lock<std::mutex> guard(lock);
	changed.wait(guard, [&] { return submitted - taken < OTHER_JOBS; });
	const long long number = submitted;
	Job& j = jobs[number % OTHER_JOBS];
	guard.unlock();

	gpu_start();
	Reserve(j, (size_t)count + count / 4);
	j.count = count;
	for (int i = 0; i < count && a != NULL; i++)
	{
		j.host_ab[2 * (size_t)i] = a[i];
		j.host_ab[2 * (size_t)i + 1] = b[i];
	}
	if (a == NULL && count > 0)
	{
		memcpy(j.host_values.data(), limbs, (size_t)count * FACTOR_LIMBS * sizeof(uint64_t));
		memcpy(j.host_sizes.data(), sizes, (size_t)count * sizeof(uint32_t));
	}
	j.host_header[0] = (uint32_t)count;
	j.host_header[1] = 0;
	if (count > 0)
	{
		const uint32_t pool = (uint32_t)count * FACTOR_MAX;
		const unsigned int blocks = (unsigned int)std::min(1024, (count + FACTOR_BLOCK / 32 - 1) / (FACTOR_BLOCK / 32));
		CUDA_CHECK(cudaMemcpyAsync(j.header, j.host_header.data(), 2 * sizeof(uint32_t), cudaMemcpyHostToDevice, j.stream));
		if (a != NULL)
		{
			CUDA_CHECK(cudaMemcpyAsync(j.ab, j.host_ab.data(), 2 * (size_t)count * sizeof(int64_t), cudaMemcpyHostToDevice, j.stream));
			k_other_norms<<<(count + 127) / 128, 128, 0, j.stream>>>(j.list, j.header, j.ab, coefficients, degree, coefficient_bits, j.values, j.sizes);
		}
		else
		{
			CUDA_CHECK(cudaMemcpyAsync(j.values, j.host_values.data(), (size_t)count * FACTOR_LIMBS * sizeof(uint64_t), cudaMemcpyHostToDevice, j.stream));
			CUDA_CHECK(cudaMemcpyAsync(j.sizes, j.host_sizes.data(), (size_t)count * sizeof(uint32_t), cudaMemcpyHostToDevice, j.stream));
		}
		k_factor_smooth<<<blocks, FACTOR_BLOCK, 0, j.stream>>>(j.list, j.header, j.values, FACTOR_LIMBS, j.sizes, NULL, j.counts, j.offsets, j.primes,
			j.header + 1, pool);
		CUDA_CHECK(cudaGetLastError());
		CUDA_CHECK(cudaMemcpyAsync(j.host_counts.data(), j.counts, (size_t)count * sizeof(int), cudaMemcpyDeviceToHost, j.stream));
		CUDA_CHECK(cudaMemcpyAsync(j.host_offsets.data(), j.offsets, (size_t)count * sizeof(uint32_t), cudaMemcpyDeviceToHost, j.stream));
		CUDA_CHECK(cudaMemcpyAsync(j.host_primes.data(), j.primes, (size_t)pool * sizeof(uint32_t), cudaMemcpyDeviceToHost, j.stream));
	}
	CUDA_CHECK(cudaEventRecord(j.done, j.stream));

	guard.lock();
	j.busy = true;
	submitted++;
	return (int)(number % (1LL << 30));
}

//The factors of the norms of a job, in the order they were submitted: for norm i, counts[i] of them (or -1 if it was not factored)
//from offsets[i] in primes, in increasing order. The jobs must be taken in the order they were submitted.
void Gpu_other_side_factoring::Results(int job, std::vector<int>& counts, std::vector<unsigned int>& offsets, std::vector<unsigned int>& primes)
{
	Job& j = jobs[(taken % (1LL << 30)) % OTHER_JOBS];
	if (job != (int)(taken % (1LL << 30)) || !j.busy)
	{
		fprintf(stderr, "Gpu_other_side_factoring: results of job %d taken out of order\n", job);
		exit(3);
	}
	gpu_start();
	CUDA_CHECK(cudaEventSynchronize(j.done));
	counts.assign(j.host_counts.data(), j.host_counts.data() + j.count);
	offsets.assign(j.host_offsets.data(), j.host_offsets.data() + j.count);
	primes.assign(j.host_primes.data(), j.host_primes.data() + (size_t)j.count * FACTOR_MAX);
	{
		std::lock_guard<std::mutex> guard(lock);
		j.busy = false;
		taken++;
	}
	changed.notify_all();
}

//The polynomial of the other side, as numbers of OTHER_WIDTH limbs in two's complement, if it is small enough: k_other_norms then
//computes the norms for a and b of up to (64 OTHER_WIDTH - 5 - coefficient bits) / degree bits.
Other_side_factoring* create_gpu_other_side_factoring(const std::vector<mpz_class>& coefficients)
{
	const int degree = (int)coefficients.size() - 1;
	if (degree < 1 || degree > OTHER_MAX_DEGREE)
	{
		return NULL;
	}
	int bits = 0;
	std::vector<uint64_t> words((size_t)(degree + 1) * OTHER_WIDTH, 0);
	for (int e = 0; e <= degree; e++)
	{
		mpz_srcptr c = coefficients[e].get_mpz_t();
		const int c_bits = (int)mpz_sizeinbase(c, 2);
		if (c_bits + degree * 32 + 4 > 64 * OTHER_WIDTH - 1)
		{
			return NULL; //too large even for a and b of 32 bits
		}
		bits = (c_bits > bits) ? c_bits : bits;
		uint64_t* x = &words[(size_t)e * OTHER_WIDTH];
		for (size_t k = 0; k < mpz_size(c); k++)
		{
			x[k] = mpz_getlimbn(c, k);
		}
		if (mpz_sgn(c) < 0)
		{
			uint64_t carry = 1;
			for (int k = 0; k < OTHER_WIDTH; k++)
			{
				x[k] = ~x[k] + carry;
				carry = (x[k] < carry);
			}
		}
	}
	return new Gpu_other_side_factoring(words, degree, bits);
}

//The GPU used, chosen by gpu_init. Each host thread has its own current device, so the threads that use the GPU set it first
//(gpu_start).
static int gpu_device = -1;
static thread_local bool gpu_device_set = false;

bool gpu_init(int device)
{
	int count = 0;
	cudaError_t error = cudaGetDeviceCount(&count);
	if (error != cudaSuccess || count == 0)
	{
		printf("This is the GPU build of NFS_factory, but no CUDA GPU was found (%s).\n", (error != cudaSuccess) ? cudaGetErrorString(error) : "no devices");
		printf("Use the CPU build (NFS_factory) on this computer, or set backend: mpir in the polynomial file.\n");
		return false;
	}
	if (device < 0 || device >= count)
	{
		printf("There is no GPU %d: this computer has %d CUDA GPU%s, numbered from 0.\n", device, count, (count > 1) ? "s" : "");
		return false;
	}
	cudaDeviceProp properties;
	CUDA_CHECK(cudaGetDeviceProperties(&properties, device));
	if (properties.major * 10 + properties.minor < 75)
	{
		printf("GPU %d (%s) has compute capability %d.%d, and the GPU build needs 7.5 or later. Use the CPU build (NFS_factory).\n", device, properties.name,
			properties.major, properties.minor);
		return false;
	}
	gpu_device = device;
	CUDA_CHECK(cudaSetDevice(device));
	//The threads that wait for the GPU sleep instead of spinning, which leaves the CPU to the threads that check the batches
	//(the waits get a little longer).
	CUDA_CHECK(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
	gpu_device_set = true;
	printf("GPU %d: %s, compute capability %d.%d, %.1f GB.\n", device, properties.name, properties.major, properties.minor, properties.totalGlobalMem / 1e9);
	return true;
}

//Before a thread uses the GPU.
static void gpu_start()
{
	if (gpu_device < 0 && !gpu_init(0))
	{
		exit(3);
	}
	if (!gpu_device_set)
	{
		CUDA_CHECK(cudaSetDevice(gpu_device));
		gpu_device_set = true;
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
			gpu.Combine_products(plan, PRODUCTS_ADJACENT, NULL, X, Y, pairs, false, parents, 0);
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

//As many trees as fit in 70 % of the free memory of the GPU (the rest is left for the second stage, the products of the setup,
//and other programs), from 1 to 5. A tree takes about 20 to 25 bytes per bit of its leaves, measured for 17 to 21 levels.
int gpu_pipeline(int levels, int max_bits_in_input, bool report)
{
	gpu_start();
	size_t free_bytes = 0, total_bytes = 0;
	if (cudaMemGetInfo(&free_bytes, &total_bytes) != cudaSuccess)
	{
		return 2;
	}
	const double tree_bytes = 26.0 * (double)((size_t)1 << (levels - 1)) * max_bits_in_input + 64e6;
	int trees = (int)(0.7 * (double)free_bytes / tree_bytes);
	if (trees < 1 && report)
	{
		printf("GPU: %.1f GB free may not be enough for a tree of about %.2f GB; treelevels in the polynomial file can make it smaller.\n", free_bytes / 1e9, tree_bytes / 1e9);
	}
	trees = (trees < 1) ? 1 : ((trees > GPU_MAX_PIPELINE) ? GPU_MAX_PIPELINE : trees);
	if (report)
	{
		printf("GPU: %.1f of %.1f GB free, about %.2f GB per tree: %d batches in progress at once.\n", free_bytes / 1e9, total_bytes / 1e9, tree_bytes / 1e9, trees);
	}
	return trees;
}

Tree_backend* create_gpu_tree_backend(int levels, int max_bits_in_input)
{
	gpu_start();
	return new Gpu_tree_backend(levels, max_bits_in_input);
}
