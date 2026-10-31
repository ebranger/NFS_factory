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

private:
	Gpu_ntt_tables tables;
	void* table_memory;
	uint64_t* scratch;
	size_t scratch_limbs;
	uint8_t* carry_scratch; //the carry functions of the limbs and chunks of Combine
	size_t carry_bytes;
	void Reserve_scratch(size_t limbs);
};
