//The integer NTT on the GPU: the transforms of ntt_core and the steps of the products of int_ntt, for batches of integers that
//share a plan. Every result is bit for bit the same as on the CPU.
//
//Integers are stored in levels: count numbers of up to stride 64-bit limbs each, number r at limbs + r stride, with size[r] limbs
//without leading zeros (limbs above the size are not read). Transforms of rows numbers are stored as in gpuecm2:
//X[(q rows + row) T + i] for prime q.
#pragma once
#include <stdint.h>
#include <stddef.h>
#include <cuda_runtime.h>
#include "int_ntt.h"

struct Gpu_level
{
	uint64_t* limbs;
	size_t stride;
	uint32_t* size;
	size_t count;
};

//The tables of the primes on the GPU, passed to the kernels by value.
struct Gpu_ntt_tables
{
	int log_max; //NTT_MAX_LOG_LENGTH
	uint32_t p[NTT_MAX_PRIMES], p_inv[NTT_MAX_PRIMES];
	uint32_t two32[NTT_MAX_PRIMES], two32_shoup[NTT_MAX_PRIMES];
	uint32_t scale[NTT_MAX_PRIMES][NTT_MAX_LOG_LENGTH + 1], scale_shoup[NTT_MAX_PRIMES][NTT_MAX_LOG_LENGTH + 1];
	uint32_t garner[NTT_MAX_PRIMES][NTT_MAX_PRIMES], garner_shoup[NTT_MAX_PRIMES][NTT_MAX_PRIMES];
	//NTT twiddles as in gpuecm2: w = omega_4096^t for t < 2048, with Shoup factors, forward and inverse, and the twist tables in
	//Montgomery form: TL = omega_Tmax^t (t < 2048), TH = omega_Tmax^(2048 h) (h < Tmax / 2048).
	const uint32_t *w, *wp, *iw, *iwp; //[prime][2048]
	const uint32_t *TL, *TH, *iTL, *iTH; //[prime][2048], [prime][Tmax / 2048]
};

#define CUDA_CHECK(call) cuda_check((call), #call, __FILE__, __LINE__)
void cuda_check(cudaError_t error, const char* call, const char* file, int line);

class Gpu_ntt
{
public:
	Gpu_ntt();
	~Gpu_ntt();

	//The limbs a level needs for the integers that Combine gives with this plan.
	static size_t Combine_limbs(const Ntt_plan& plan);

	//The transforms of the numbers of a level (each less than 2^L): primes rows of count transforms.
	void Forward(const Ntt_plan& plan, const Gpu_level& in, uint32_t* X, cudaStream_t stream);

	//Forward and inverse transforms of rows rows, in place, as ntt_forward and ntt_inverse.
	void Transform(const Ntt_plan& plan, uint32_t* X, size_t rows, cudaStream_t stream);
	void Inverse(const Ntt_plan& plan, uint32_t* X, size_t rows, cudaStream_t stream);

	//Pointwise products in Montgomery form, as ntt_pointwise. X = X * Y for rows rows.
	void Pointwise(const Ntt_plan& plan, uint32_t* X, const uint32_t* Y, size_t rows, cudaStream_t stream);
	//Y[j] = X[2j] * X[2j + 1] for j < count: the products of the pairs of the product tree.
	void Pointwise_adjacent(const Ntt_plan& plan, const uint32_t* X, uint32_t* Y, size_t count, cudaStream_t stream);
	//Z[2j] = B[j] * C[2j + 1] and Z[2j + 1] = B[j] * C[2j] for j < count: a parent times the sibling, in the remainder tree.
	void Pointwise_swap(const Ntt_plan& plan, const uint32_t* B, const uint32_t* C, uint32_t* Z, size_t count, cudaStream_t stream);

	//The integers from rows pointwise products, as Ntt_product::Inverse: the inverse transforms (X is overwritten), the CRT and
	//the carries, into out, whose stride must be at least Combine_limbs. Exact products, or products modulo 2^L - 1 from 0 to
	//2^L - 2 if cyclic. Uses scratch space of rows T 3 limbs, which the class keeps.
	void Combine(const Ntt_plan& plan, uint32_t* X, size_t rows, bool cyclic, Gpu_level& out, cudaStream_t stream);

	//out[r] = bits shift[r] up to precision[r] of in[r], minus one (modulo 2^(precision - shift)), as Ntt_product::Extract_middle,
	//or 0 if y_size[r / y_div] is 0 (a middle product of y = 0).
	void Extract_middle(const Gpu_level& in, const uint32_t* shift, const uint32_t* precision, const uint32_t* y_size, int y_div, Gpu_level& out, cudaStream_t stream);

	//Exact products with schoolbook multiplication, for small numbers: out[j] = x[2j] * x[2j + 1]. out's stride must be at least
	//twice x's.
	void School_products(const Gpu_level& x, Gpu_level& out, cudaStream_t stream);

	//Exact truncated products with schoolbook multiplication: out[2j] = ((y[j] * x[2j + 1]) mod 2^precision[2j]) >> shift[2j], and
	//out[2j + 1] the same with x[2j]. out's stride must be enough for the result.
	void School_middle(const Gpu_level& y, const Gpu_level& x, const uint32_t* shift, const uint32_t* precision, Gpu_level& out, cudaStream_t stream);

	//The size in bits of each number of a level (1 for 0, as mpz_sizeinbase).
	void Bits(const Gpu_level& in, uint32_t* bits, cudaStream_t stream);

	//Single numbers in device memory, as a number of limbs (zero above the value). out may not overlap in, except as noted.
	//out = floor(in / 2^position) modulo 2^count (no limit if count is 0), in out_limbs limbs. position < 0 shifts to the left.
	//With complement, of ~in instead, with ones above in_limbs: ~in modulo 2^count is 2^count - 1 - in.
	void Shift(uint64_t* out, size_t out_limbs, const uint64_t* in, size_t in_limbs, long long position, uint64_t count, bool complement, cudaStream_t stream);
	//out = a + b, or a - b, modulo 2^(64 n). a has a_limbs limbs (NULL for 0). b has b_limbs limbs at limb b_offset, or if b is
	//NULL, b is b_value at limb b_offset. The carry (borrow) out goes to *carry on the device, if not NULL. out can be a, and b if
	//b_offset is 0.
	void Add(uint64_t* out, size_t n, const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs, size_t b_offset, uint64_t b_value, bool subtract,
		uint32_t* carry, cudaStream_t stream);
	//*limbs = the number of limbs of in without leading zeros, and *bits its size in bits (1 for 0, as mpz_sizeinbase) if not NULL.
	void Size(const uint64_t* in, size_t n, uint32_t* limbs, uint32_t* bits, cudaStream_t stream);
	void Set_word(uint32_t* p, uint32_t value, cudaStream_t stream);

	//The cofactors of the leaves x of a scaled remainder tree, one thread per leaf, as Cpu_tree_backend::Leaf_cofactors: leaf i,
	//N, has the fraction of its parent at y[i / 2] (scaled by 2^(y_bits[i / 2] + guard), with y_bits the bits of the parents in the
	//product tree) and its sibling at x[i ^ 1]. out[i] = N / gcd(N, (P mod N)^(2^squarings)). The leaves have at most max_limbs
	//limbs, and out's stride must be at least Leaf_kernel_limbs(max_limbs). False, with nothing done, if the leaves are too large
	//(Leaf_kernel_limbs is 0) or squarings more than 2.
	bool Leaf_cofactors(const Gpu_level& y, const uint32_t* y_bits, const Gpu_level& x, const uint32_t* x_bits, int max_limbs, uint32_t guard, int squarings, Gpu_level& out,
		cudaStream_t stream);
	static int Leaf_kernel_limbs(int max_limbs);

private:
	Gpu_ntt_tables tables;
	void* table_memory;
	uint64_t* scratch;
	size_t scratch_limbs;
	uint8_t* carry_scratch; //the carry functions of the limbs and chunks of Combine
	size_t carry_bytes;
	unsigned int* row_done; //the blocks done in each row, for the carries
	size_t row_done_count;
	unsigned int* size_found; //the top found by the blocks of k_size, and the blocks done, 0 between calls
	size_t combine_coefficients; //that the combine keeps in shared memory, from the shared memory of the GPU
	void Reserve_scratch(size_t limbs);
	void Reserve_carry(size_t bytes);
	void Reserve_rows(size_t rows);
};
