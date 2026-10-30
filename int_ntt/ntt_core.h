//Number-theoretic transforms modulo 31-bit primes p = c 2^k + 1, with the algorithm and output order of the NTT of gpuecm2
//(src/gpu/ntt.cu), so that the CPU and GPU versions give exactly the same transforms.
//
//The forward transform of length T = 2^logT is a decimation-in-frequency NTT with omega_T = g^((p - 1) / T), where g is the
//least quadratic nonresidue from 3 upwards. Its output is in bit-reversed order: position j holds the sum of x_n omega_T^(n k)
//with k = bitreverse(j). Transforms of up to 2^NTT_ROW_LOG_LENGTH points are radix-2 transforms of the whole row. Longer ones
//use the four-step method T = R C with R = 2^max(logT - 12, 3): R-point transforms down the columns, a twist of row i and
//column l by omega_T^(l bitreverse(i)), then C-point transforms along the rows. That too gives the bit-reversed order.
//The inverse transform undoes the forward transform exactly, but times T (decimation in time, without the 1/T).
#pragma once
#include <stdint.h>
#include <stddef.h>
#include <vector>

#define NTT_MIN_LOG_LENGTH 3
#define NTT_MAX_LOG_LENGTH 24
#define NTT_ROW_LOG_LENGTH 12 //The longest transform done as a single row (gpuecm2's tile of 4096 points).
#define NTT_COLUMN_BLOCK 16 //Columns of a four-step transform done together on the CPU, one cache line per row.

//A prime p with 2^30 < p < 2^31 and 2^NTT_MAX_LOG_LENGTH dividing p - 1 (Init assumes both), and its tables.
struct Ntt_prime
{
	uint32_t p;
	uint32_t p_inv; //-1 / p mod 2^32
	uint32_t generator; //the least quadratic nonresidue from 3 upwards
	uint32_t two32, two32_shoup; //2^32 mod p

	//The twiddles of the radix-2 levels: at s + j, omega_(2s)^j for s = 1, 2, 4, ..., 2^(NTT_ROW_LOG_LENGTH - 1) and j < s,
	//with their Shoup factors. The same with the inverse roots.
	std::vector<uint32_t> w, w_shoup, iw, iw_shoup;

	//The twist factors of the four-step transforms in Montgomery form, as in gpuecm2: twist_low[t] = omega_Tmax^t for t < 2048,
	//and twist_high[h] = omega_Tmax^(2048 h) for h < Tmax / 2048, with Tmax = 2^NTT_MAX_LOG_LENGTH. The same with the inverse roots.
	std::vector<uint32_t> twist_low, twist_high, itwist_low, itwist_high;

	void Init(uint32_t prime);
};

//The scratch space that the transforms of length 2^logT need, in words.
size_t ntt_scratch_size(int logT);

//Forward transform of the 2^logT residues x, in place, with output in bit-reversed order.
void ntt_forward(uint32_t* x, int logT, const Ntt_prime& prime, uint32_t* scratch);

//Inverse of ntt_forward, times 2^logT, in place.
void ntt_inverse(uint32_t* x, int logT, const Ntt_prime& prime, uint32_t* scratch);

//x = x * y / 2^32 mod p (Montgomery) for n residues, as the pointwise products of gpuecm2.
void ntt_pointwise(uint32_t* x, const uint32_t* y, size_t n, const Ntt_prime& prime);
