//Word-level arithmetic modulo the 31-bit NTT primes, 2^30 < p < 2^31. These are the portable versions of the functions of
//arith32.cuh in gpuecm2 (by the same author), and the CUDA build of int_ntt is to share them.
//Every result is fully reduced, from 0 to p - 1, so any correct order of the same operations gives exactly the same numbers.
#pragma once
#include <stdint.h>

#if defined(__CUDACC__)
#define INT_NTT_FUNC __host__ __device__ __forceinline__
#else
#define INT_NTT_FUNC static inline
#endif

//The high 32 bits of a * b.
INT_NTT_FUNC uint32_t mul_hi32(uint32_t a, uint32_t b)
{
#if defined(__CUDA_ARCH__)
	return __umulhi(a, b);
#else
	return (uint32_t)(((uint64_t)a * b) >> 32);
#endif
}

//a + b mod p, for a and b less than p.
INT_NTT_FUNC uint32_t add_mod(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t s = a + b; //less than 2^32, as p < 2^31
	return (s >= p) ? s - p : s;
}

//a - b mod p, for a and b less than p.
INT_NTT_FUNC uint32_t sub_mod(uint32_t a, uint32_t b, uint32_t p)
{
	return (a >= b) ? a - b : a + p - b;
}

//b * w mod p for any 32-bit b and w less than p, with w_shoup = floor(w 2^32 / p) precomputed (Shoup).
INT_NTT_FUNC uint32_t mul_shoup(uint32_t b, uint32_t w, uint32_t w_shoup, uint32_t p)
{
	uint32_t q = mul_hi32(w_shoup, b);
	uint32_t r = w * b - q * p; //from 0 to 2p - 1
	return (r >= p) ? r - p : r;
}

//a * b / 2^32 mod p (Montgomery), for a * b less than 2^32 p, with p_inv = -1 / p mod 2^32.
INT_NTT_FUNC uint32_t mul_mont(uint32_t a, uint32_t b, uint32_t p, uint32_t p_inv)
{
	uint32_t lo = a * b, hi = mul_hi32(a, b), m = lo * p_inv;
	uint32_t r = hi + mul_hi32(m, p) + (lo != 0); //less than 2p
	return (r >= p) ? r - p : r;
}

//The Shoup factor floor(w 2^32 / p) of w, for w less than p.
INT_NTT_FUNC uint32_t shoup_factor(uint32_t w, uint32_t p)
{
	return (uint32_t)(((uint64_t)w << 32) / p);
}
