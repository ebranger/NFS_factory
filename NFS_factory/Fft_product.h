#pragma once
#include "mpirxx.h"
#include <vector>

//Products modulo 2^L - 1 with the FFT of MPIR, for the scaled remainder tree. They use the cache friendly matrix Fourier
//algorithm of MPIR's large products. The products modulo B^n + 1 that MPIR exports use a plain radix 2 FFT, which is
//up to twice as slow for large numbers.
class Cyclic_product
{
public:
	//result = ((y * x) mod 2^precision) >> shift, or one less (modulo 2^(precision - shift)).
	//y must be less than 2^precision, and x less than 2^(shift + 1).
	void Middle_product(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t shift, mp_bitcnt_t precision);

	//Two middle products with the same y, which is transformed only once: result0 from y and x0 with shift0, and result1 from y
	//and x1 with shift1.
	void Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision);

	//result = y * x modulo 2^L - 1, from 0 to 2^L - 2, with L at least min_L. Both numbers must be less than 2^L. Returns L.
	mp_bitcnt_t Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L);

	//An approximation of B^(2n) / a, where a has n limbs and its top bit set, with B = 2^64. The result is at most
	//INVERSE_MAX_ERROR less than floor((B^(2n) - 1) / a), and never more.
	void Inverse(mpz_ptr result, mpz_srcptr a);

private:
	void Newton_inverse(mpz_ptr result, mpz_srcptr a, mp_size_t a_limbs);
	void Setup(mp_size_t y_limbs, mp_size_t x_limbs, mp_bitcnt_t min_L);
	void Transform_columns(mp_ptr* c, mp_srcptr a, mp_size_t a_limbs);
	void Transform_rows(mp_ptr* c);
	void Multiply_rows(mp_ptr* x, mp_ptr* y);
	void Fold(mp_ptr* c);
	mp_bitcnt_t Fold_product(mp_srcptr y, mp_size_t y_limbs, mp_srcptr x, mp_size_t x_limbs, mp_bitcnt_t min_L);
	void Extract_middle(mpz_ptr result, mp_bitcnt_t shift, mp_bitcnt_t precision);

	//The transform: 4n = 2^(depth + 2) coefficients of bits1 bits, computed modulo 2^(n w) + 1 in limbs + 1 limbs, for products
	//modulo 2^L - 1 with L = 4n bits1. The matrix Fourier algorithm splits the transform into rows of n1 coefficients.
	int depth;
	mp_size_t n, w, limbs, n1;
	mp_bitcnt_t bits1, L;

	std::vector<mp_limb_t> coefficients; //the coefficients of both transforms, and scratch space for the FFT.
	std::vector<mp_ptr> ii, jj; //pointers to the coefficients.
	mp_ptr t1, t2, s1, tt; //scratch space.
	std::vector<mp_limb_t> product; //the product modulo 2^L - 1, folded to less than 2^(L + 1).
	mp_size_t product_limbs;
	mpz_class high_part;
};

#define INVERSE_MAX_ERROR (1 << 16)
