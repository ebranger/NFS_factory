//Products modulo 2^L - 1 with the matrix Fourier algorithm FFT of MPIR.
//The numbers are split into 4n coefficients of bits1 bits, and a transform of length 4n computes their cyclic convolution with
//coefficients modulo 2^(n w) + 1, where 2^w is a 2n-th root of unity (and the square root of 2 trick gives a 4n-th root).
//A cyclic convolution of length 4n is a product modulo 2^L - 1 with L = 4n bits1, once the coefficients are added together
//and the part above L is folded back onto the bottom.
//The matrix Fourier algorithm does the transform as FFTs of length 2n / n1 on the columns of a matrix with rows of n1 coefficients,
//and then FFTs of length n1 on the rows.

#include "Fft_product.h"
#include <stdlib.h>

#define FFT_MAX_COEFFICIENT_LIMBS 256 //Larger coefficients would need the sizes that MPIR's pointwise products modulo B^n + 1 accept.
#define FFT_MIN_DEPTH 4
#define INVERSE_EXACT_LIMBS 2000 //Inverses of up to this many limbs are computed exactly by MPIR.
#define INVERSE_MARGIN (1 << 12) //Subtracted from the Newton inverse, which can be a few units too large.

//MPIR functions for the FFT. They are declared in gmp-impl.h, which is not part of the installed headers.
extern "C"
{
	mp_size_t __mpir_fft_split_bits(mp_ptr* poly, mp_srcptr limbs, mp_size_t total_limbs, mp_bitcnt_t bits, mp_size_t output_limbs);
	void __mpir_fft_combine_bits(mp_ptr res, const mp_ptr* poly, long length, mp_bitcnt_t bits, mp_size_t output_limbs, mp_size_t total_limbs);
	void __mpir_fft_mfa_trunc_sqrt2_outer(mp_ptr* ii, mp_size_t n, mp_bitcnt_t w, mp_ptr* t1, mp_ptr* t2, mp_ptr* temp, mp_size_t n1, mp_size_t trunc);
	void __mpir_ifft_mfa_trunc_sqrt2_outer(mp_ptr* ii, mp_size_t n, mp_bitcnt_t w, mp_ptr* t1, mp_ptr* t2, mp_ptr* temp, mp_size_t n1, mp_size_t trunc);
	void __mpir_fft_mfa_trunc_sqrt2_inner(mp_ptr* ii, mp_ptr* jj, mp_size_t n, mp_bitcnt_t w, mp_ptr* t1, mp_ptr* t2, mp_ptr* temp, mp_size_t n1, mp_size_t trunc, mp_ptr tt);
	void __mpir_fft_radix2(mp_ptr* ii, mp_size_t n, mp_bitcnt_t w, mp_ptr* t1, mp_ptr* t2);
	void __mpir_ifft_radix2(mp_ptr* ii, mp_size_t n, mp_bitcnt_t w, mp_ptr* t1, mp_ptr* t2);
	void __gmpn_normmod_2expp1(mp_ptr t, mp_size_t l);
	int __gmpn_mulmod_Bexpp1(mp_ptr r, mp_srcptr i1, mp_srcptr i2, mp_size_t limbs, mp_ptr tt);
}

//The fastest depth (n = 2^depth) for a product modulo 2^L - 1 of L limbs is the first depth whose limit is above L.
//Tuned on an i7-14700F.
static const mp_size_t depth_limits[] = { 1400, 2200, 7000, 14000, 60000, 110000, 230000, 380000, 800000 };

//The transform parameters for a product of numbers of y_limbs and x_limbs limbs modulo 2^L - 1, with L at least min_L.
//Returns false if there are none at this depth.
static bool choose_parameters(int depth, mp_size_t y_limbs, mp_size_t x_limbs, mp_bitcnt_t min_L, mp_size_t& w, mp_bitcnt_t& bits1)
{
	mp_size_t n = (mp_size_t)1 << depth;
	mp_size_t w_step = (n >= GMP_LIMB_BITS) ? 1 : GMP_LIMB_BITS / n; //n w must be a whole number of limbs.
	mp_bitcnt_t min_bits1 = (min_L + 4 * n - 1) / (4 * n);
	w = (2 * min_bits1 + depth + 1 + n - 1) / n;
	w = (w + w_step - 1) / w_step * w_step;

	for (;; w += w_step)
	{
		if (n * w / GMP_LIMB_BITS > FFT_MAX_COEFFICIENT_LIMBS)
		{
			return false;
		}

		//Each coefficient of the convolution is a sum of products of two coefficients, one for each coefficient of the shorter
		//number, and the sum must fit in n w bits: up to 2n products with bits1 = (n w - depth - 1) / 2, and up to 4n with one bit less.
		bits1 = (n * w - (depth + 1)) / 2;
		mp_size_t j1 = (y_limbs * GMP_LIMB_BITS - 1) / bits1 + 1;
		mp_size_t j2 = (x_limbs * GMP_LIMB_BITS - 1) / bits1 + 1;

		if (j1 > 2 * n && j2 > 2 * n)
		{
			bits1 = (n * w - (depth + 2)) / 2;
			j1 = (y_limbs * GMP_LIMB_BITS - 1) / bits1 + 1;
			j2 = (x_limbs * GMP_LIMB_BITS - 1) / bits1 + 1;
		}

		if (4 * n * bits1 >= min_L && j1 <= 4 * n && j2 <= 4 * n)
		{
			return true;
		}
	}
}

//Choose the transform for products of numbers of up to y_limbs and x_limbs limbs modulo 2^L - 1, with L at least min_L,
//and set up the coefficients.
void Cyclic_product::Setup(mp_size_t y_limbs, mp_size_t x_limbs, mp_bitcnt_t min_L)
{
	depth = FFT_MIN_DEPTH;
	while (depth - FFT_MIN_DEPTH < (int)(sizeof(depth_limits) / sizeof(depth_limits[0])) && (mp_size_t)(min_L / GMP_LIMB_BITS) >= depth_limits[depth - FFT_MIN_DEPTH])
	{
		depth++;
	}

	while (!choose_parameters(depth, y_limbs, x_limbs, min_L, w, bits1))
	{
		depth++;
	}

	n = (mp_size_t)1 << depth;
	limbs = (n * w) / GMP_LIMB_BITS;
	n1 = (mp_size_t)1 << (depth / 2);
	L = 4 * n * bits1;
	mp_size_t size = limbs + 1;

	//Space for the coefficients of both numbers, and for the temporaries of the FFT: three coefficients, and two for the pointwise products.
	if (coefficients.size() < (size_t)(8 * n * size + 5 * size))
	{
		coefficients.resize(8 * n * size + 5 * size);
	}
	ii.resize(4 * n);
	jj.resize(4 * n);
	mp_ptr ptr = coefficients.data();
	for (mp_size_t i = 0; i < 4 * n; i++, ptr += size)
	{
		ii[i] = ptr;
	}
	for (mp_size_t i = 0; i < 4 * n; i++, ptr += size)
	{
		jj[i] = ptr;
	}
	t1 = ptr;
	t2 = t1 + size;
	s1 = t2 + size;
	tt = s1 + size;
}

//Split a number into the coefficients c, and do the FFTs on the columns.
void Cyclic_product::Transform_columns(mp_ptr* c, mp_srcptr a, mp_size_t a_limbs)
{
	mp_size_t j = __mpir_fft_split_bits(c, a, a_limbs, bits1, limbs);
	for (; j < 4 * n; j++)
	{
		mpn_zero(c[j], limbs + 1);
	}
	__mpir_fft_mfa_trunc_sqrt2_outer(c, n, w, &t1, &t2, &s1, n1, 4 * n);
}

//Complete the transform of c, after Transform_columns, with the FFTs on the rows.
void Cyclic_product::Transform_rows(mp_ptr* c)
{
	mp_size_t rows = 4 * n / n1;
	for (mp_size_t r = 0; r < rows; r++)
	{
		__mpir_fft_radix2(c + r * n1, n1 / 2, w * (2 * n / n1), &t1, &t2);
		for (mp_size_t j = 0; j < n1; j++)
		{
			__gmpn_normmod_2expp1(c[r * n1 + j], limbs);
		}
	}
}

//Multiply x, after Transform_columns, with the fully transformed y, and transform the product back. One row at a time, this does
//the FFT of the row of x, the pointwise products, and the inverse FFT of the row, as MPIR does for its products.
void Cyclic_product::Multiply_rows(mp_ptr* x, mp_ptr* y)
{
	mp_size_t rows = 4 * n / n1;
	for (mp_size_t r = 0; r < rows; r++)
	{
		__mpir_fft_radix2(x + r * n1, n1 / 2, w * (2 * n / n1), &t1, &t2);
		for (mp_size_t j = 0; j < n1; j++)
		{
			mp_size_t t = r * n1 + j;
			__gmpn_normmod_2expp1(x[t], limbs);
			__gmpn_mulmod_Bexpp1(x[t], y[t], x[t], limbs, tt);
		}
		__mpir_ifft_radix2(x + r * n1, n1 / 2, w * (2 * n / n1), &t1, &t2);
	}
	__mpir_ifft_mfa_trunc_sqrt2_outer(x, n, w, &t1, &t2, &s1, n1, 4 * n);
}

//Add the coefficients c of a product together into product, and fold the part from L upwards back onto the bottom. The result is
//the product modulo 2^L - 1, less than 2^(L + 1).
void Cyclic_product::Fold(mp_ptr* c)
{
	//The last coefficients reach up to n w bits past L.
	mp_size_t L_limbs = (mp_size_t)(L / GMP_LIMB_BITS);
	unsigned int L_bit = (unsigned int)(L % GMP_LIMB_BITS);
	product_limbs = L_limbs + limbs + 2;
	mp_size_t high_limbs = product_limbs - L_limbs;
	product.assign(product_limbs + high_limbs, 0);
	mp_ptr p = product.data();
	__mpir_fft_combine_bits(p, c, (long)(4 * n), bits1, limbs, product_limbs);

	mp_ptr high = p + product_limbs;
	if (L_bit != 0)
	{
		mpn_rshift(high, p + L_limbs, high_limbs, L_bit);
		p[L_limbs] &= (((mp_limb_t)1) << L_bit) - 1;
		mpn_zero(p + L_limbs + 1, high_limbs - 1);
	}
	else
	{
		mpn_copyi(high, p + L_limbs, high_limbs);
		mpn_zero(p + L_limbs, high_limbs);
	}
	while (high_limbs > 0 && high[high_limbs - 1] == 0)
	{
		high_limbs--;
	}
	if (high_limbs > 0)
	{
		mpn_add(p, p, product_limbs, high, high_limbs);
	}
}

//Compute y * x modulo 2^L - 1, with L at least min_L, both numbers less than 2^L, and neither zero. The result is put in
//product, less than 2^(L + 1), and L is returned.
mp_bitcnt_t Cyclic_product::Fold_product(mp_srcptr y, mp_size_t y_limbs, mp_srcptr x, mp_size_t x_limbs, mp_bitcnt_t min_L)
{
	Setup(y_limbs, x_limbs, min_L);
	Transform_columns(ii.data(), y, y_limbs);
	Transform_columns(jj.data(), x, x_limbs);
	__mpir_fft_mfa_trunc_sqrt2_inner(ii.data(), jj.data(), n, w, &t1, &t2, &s1, n1, 4 * n, tt);
	__mpir_ifft_mfa_trunc_sqrt2_outer(ii.data(), n, w, &t1, &t2, &s1, n1, 4 * n);
	Fold(ii.data());
	return L;
}

//From the folded product, result = bits from shift up to precision, minus one.
void Cyclic_product::Extract_middle(mpz_ptr result, mp_bitcnt_t shift, mp_bitcnt_t precision)
{
	mp_ptr p = product.data();
	mp_bitcnt_t out_bits = precision - shift;
	mp_size_t out_limbs = (mp_size_t)((out_bits + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS);
	mp_size_t first_limb = (mp_size_t)(shift / GMP_LIMB_BITS);
	mp_size_t source_limbs = (mp_size_t)((precision + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS) - first_limb;

	if (result->_mp_alloc < source_limbs)
	{
		mpz_realloc2(result, source_limbs * GMP_LIMB_BITS);
	}

	if (shift % GMP_LIMB_BITS != 0)
	{
		mpn_rshift(result->_mp_d, p + first_limb, source_limbs, (unsigned int)(shift % GMP_LIMB_BITS));
	}
	else
	{
		mpn_copyi(result->_mp_d, p + first_limb, source_limbs);
	}
	mpn_sub_1(result->_mp_d, result->_mp_d, out_limbs, 1);

	if (out_bits % GMP_LIMB_BITS != 0)
	{
		result->_mp_d[out_limbs - 1] &= (((mp_limb_t)1) << (out_bits % GMP_LIMB_BITS)) - 1;
	}

	while (out_limbs > 0 && result->_mp_d[out_limbs - 1] == 0)
	{
		out_limbs--;
	}
	result->_mp_size = (int)out_limbs;
}

//result = ((y * x) mod 2^precision) >> shift, or one less (modulo 2^(precision - shift)).
//The product is computed modulo 2^L - 1, with L at least 64 bits more than precision + bits of x - shift. The part of the product
//above 2^L is then less than 2^(shift - 64). It wraps around and is added to the low bits, below the bits that are needed, which
//can make the result one too large. One is subtracted, so that the result is correct or one too small, as with the truncated products
//of the scaled remainder tree, whose guard bits absorb such errors.
void Cyclic_product::Middle_product(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t shift, mp_bitcnt_t precision)
{
	mp_size_t y_limbs = abs(y->_mp_size);
	mp_size_t x_limbs = abs(x->_mp_size);

	if (y_limbs == 0 || x_limbs == 0)
	{
		mpz_set_ui(result, 0);
		return;
	}

	Fold_product(y->_mp_d, y_limbs, x->_mp_d, x_limbs, precision + mpz_sizeinbase(x, 2) - shift + GMP_LIMB_BITS);
	Extract_middle(result, shift, precision);
}

//Two middle products with the same y, as Middle_product, but with the transform of y computed only once.
void Cyclic_product::Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision)
{
	mp_size_t y_limbs = abs(y->_mp_size);
	mp_size_t x0_limbs = abs(x0->_mp_size);
	mp_size_t x1_limbs = abs(x1->_mp_size);

	if (y_limbs == 0 || x0_limbs == 0 || x1_limbs == 0)
	{
		Middle_product(result0, y, x0, shift0, precision);
		Middle_product(result1, y, x1, shift1, precision);
		return;
	}

	mp_bitcnt_t min_L0 = precision + mpz_sizeinbase(x0, 2) - shift0 + GMP_LIMB_BITS;
	mp_bitcnt_t min_L1 = precision + mpz_sizeinbase(x1, 2) - shift1 + GMP_LIMB_BITS;
	Setup(y_limbs, (x0_limbs > x1_limbs) ? x0_limbs : x1_limbs, (min_L0 > min_L1) ? min_L0 : min_L1);

	Transform_columns(ii.data(), y->_mp_d, y_limbs);
	Transform_rows(ii.data());

	Transform_columns(jj.data(), x0->_mp_d, x0_limbs);
	Multiply_rows(jj.data(), ii.data());
	Fold(jj.data());
	Extract_middle(result0, shift0, precision);

	Transform_columns(jj.data(), x1->_mp_d, x1_limbs);
	Multiply_rows(jj.data(), ii.data());
	Fold(jj.data());
	Extract_middle(result1, shift1, precision);
}

//result = y * x modulo 2^L - 1, from 0 to 2^L - 2, with L at least min_L. Both numbers must be less than 2^L. Returns L.
mp_bitcnt_t Cyclic_product::Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L)
{
	mp_size_t y_limbs = abs(y->_mp_size);
	mp_size_t x_limbs = abs(x->_mp_size);

	if (y_limbs == 0 || x_limbs == 0)
	{
		mpz_set_ui(result, 0);
		return min_L;
	}

	Fold_product(y->_mp_d, y_limbs, x->_mp_d, x_limbs, min_L);

	mp_size_t size = product_limbs;
	while (size > 0 && product[size - 1] == 0)
	{
		size--;
	}
	if (result->_mp_alloc < size + 1)
	{
		mpz_realloc2(result, (size + 1) * GMP_LIMB_BITS);
	}
	mpn_copyi(result->_mp_d, product.data(), size);
	result->_mp_size = (int)size;

	//Fold until the result is less than 2^L, and replace 2^L - 1 by 0.
	while (mpz_sizeinbase(result, 2) > L)
	{
		mpz_tdiv_q_2exp(high_part.get_mpz_t(), result, L);
		mpz_tdiv_r_2exp(result, result, L);
		mpz_add(result, result, high_part.get_mpz_t());
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
	return L;
}

//An approximation of B^(2n) / a, where a has n limbs and its top bit set, as B^n plus n limbs. It is within a few units of
//floor((B^(2n) - 1) / a), but can be on either side of it.
//One Newton step from the inverse v_h of the top h limbs of a, about half of them: with D = a v_h - B^(n + h), which is less
//than a few B^n in size, the inverse is v_h B^(n - h) - v_h D / B^(2h). Only D modulo 2^L - 1 with L a little more than n limbs
//is needed, and the error of v_h is squared, so the result is as accurate as the truncations of the last step allow.
void Cyclic_product::Newton_inverse(mpz_ptr result, mpz_srcptr a, mp_size_t a_limbs)
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
	d *= v_h;
	mpz_fdiv_q_2exp(d.get_mpz_t(), d.get_mpz_t(), GMP_LIMB_BITS * (h + 1));

	mpz_mul_2exp(result, v_h.get_mpz_t(), GMP_LIMB_BITS * (a_limbs - h));
	mpz_sub(result, result, d.get_mpz_t());
}

//An approximation of B^(2n) / a, where a has n limbs and its top bit set, at most INVERSE_MAX_ERROR less than
//floor((B^(2n) - 1) / a), and never more.
void Cyclic_product::Inverse(mpz_ptr result, mpz_srcptr a)
{
	mp_size_t a_limbs = abs(a->_mp_size);
	Newton_inverse(result, a, a_limbs);
	if (a_limbs > INVERSE_EXACT_LIMBS)
	{
		mpz_sub_ui(result, result, INVERSE_MARGIN);
	}
}
