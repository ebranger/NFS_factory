//Integer products with the number-theoretic transforms of ntt_core.h: exact products, products modulo 2^L - 1, and the middle
//products of the scaled remainder tree. This is the CPU reference of the GPU products: both give exactly the same integers.
//
//An integer is split into coefficients of b bits (b <= 64), and the cyclic convolution of length T = 2^logT of the coefficients
//of two integers is computed modulo P primes (3 to 5). Each coefficient of the convolution is less than n 2^(2b), with n the
//number of coefficients of the shorter integer, so it is recovered exactly by the CRT when the product of the primes is larger.
//The coefficients are then added together at their bit positions. Without wraparound (T at least the number of coefficients
//of the product) that is the exact product; with wraparound it is the product modulo 2^L - 1 with L = T b, once the part from
//bit L upwards is added back onto the bottom (end-around carry).
#pragma once
#include "mpir.h"
#include "ntt_core.h"
#include <vector>

#define NTT_MIN_PRIMES 3
#define NTT_MAX_PRIMES 5
#define NTT_MAX_COEFFICIENT_BITS 64

//The transform of a product: length 2^log_length, coefficients of coefficient_bits bits, and the number of primes.
struct Ntt_plan
{
	int primes;
	int log_length;
	int coefficient_bits;

	size_t Length() const { return (size_t)1 << log_length; }
	mp_bitcnt_t L() const { return (mp_bitcnt_t)Length() * coefficient_bits; }
};

//The primes and the constants of the CRT, built on first use.
struct Ntt_tables
{
	Ntt_prime prime[NTT_MAX_PRIMES];
	int modulus_bits[NTT_MAX_PRIMES + 1]; //floor(log2) of the product of the first P primes
	uint32_t garner[NTT_MAX_PRIMES][NTT_MAX_PRIMES], garner_shoup[NTT_MAX_PRIMES][NTT_MAX_PRIMES]; //[j][i]: 1 / p_j mod p_i, j < i
	uint32_t scale[NTT_MAX_PRIMES][NTT_MAX_LOG_LENGTH + 1], scale_shoup[NTT_MAX_PRIMES][NTT_MAX_LOG_LENGTH + 1]; //2^32 / 2^logT mod p
};
const Ntt_tables& ntt_tables();

//The largest coefficient size for which convolutions with up to n terms per coefficient are exact with this many primes.
int ntt_max_coefficient_bits(int primes, size_t n);

//The cheapest plan for an exact product of numbers of a_bits and b_bits bits, both nonzero. False if there is none.
bool ntt_plan_product(Ntt_plan& plan, mp_bitcnt_t a_bits, mp_bitcnt_t b_bits);

//The cheapest plan for a product modulo 2^L - 1 with L at least min_L, of numbers of up to a_bits and b_bits bits.
//L is also at least a_bits and b_bits. False if there is none.
bool ntt_plan_cyclic(Ntt_plan& plan, mp_bitcnt_t min_L, mp_bitcnt_t a_bits, mp_bitcnt_t b_bits);

//Whether a plan gives exact results for numbers of a_bits and b_bits bits, as an exact product or modulo 2^L - 1.
bool ntt_plan_valid(const Ntt_plan& plan, mp_bitcnt_t a_bits, mp_bitcnt_t b_bits, bool cyclic);

class Ntt_product
{
public:
	//result = a * b.
	void Multiply(mpz_ptr result, mpz_srcptr a, mpz_srcptr b);

	//result = y * x modulo 2^L - 1, from 0 to 2^L - 2, with L at least min_L. Both numbers must be less than 2^L. Returns L.
	mp_bitcnt_t Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L);

	//result = ((y * x) mod 2^precision) >> shift, or one less (modulo 2^(precision - shift)), as Cyclic_product::Middle_product.
	//y must be less than 2^precision, and x less than 2^(shift + 1).
	void Middle_product(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t shift, mp_bitcnt_t precision);

	//Two middle products with the same y, which is transformed only once: result0 from y and x0 with shift0, and result1 from y
	//and x1 with shift1.
	void Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision);

	//The steps of a product, for products with a given plan and for transforms that are used more than once.
	//The transform of a number of a_limbs limbs, less than 2^L: primes rows of 2^log_length residues.
	void Forward(std::vector<uint32_t>& transform, const Ntt_plan& plan, mp_srcptr a, mp_size_t a_limbs);

	//x = x * y pointwise (in Montgomery form, as the inverse expects).
	void Pointwise(std::vector<uint32_t>& x, const std::vector<uint32_t>& y, const Ntt_plan& plan);

	//The integer from a pointwise product: the inverse transforms, the CRT and the carries. The transform is overwritten.
	//If cyclic is false, the result is the convolution as an integer, which is the exact product when the plan has no wraparound.
	//If cyclic is true, it is the product modulo 2^L - 1, from 0 to 2^L - 2. Returns the number of limbs, without leading zeros.
	mp_size_t Inverse(std::vector<mp_limb_t>& result, std::vector<uint32_t>& transform, const Ntt_plan& plan, bool cyclic);

	//Forward of a nonnegative number less than 2^L.
	void Forward(std::vector<uint32_t>& transform, const Ntt_plan& plan, mpz_srcptr a);

	//result = the exact product from a pointwise product, with a plan that has no wraparound for it. The transform is overwritten.
	void Product_from(mpz_ptr result, std::vector<uint32_t>& transform, const Ntt_plan& plan);

	//result = ((y * x) mod 2^precision) >> shift, or one less, from the pointwise product of y and x, as Middle_product. The plan
	//must have L at least max(precision + bits of x - shift, precision) + 64, and y must not be 0. The transform is overwritten.
	void Middle_from(mpz_ptr result, std::vector<uint32_t>& transform, const Ntt_plan& plan, mp_bitcnt_t shift, mp_bitcnt_t precision);

private:
	void Extract_middle(mpz_ptr result, mp_bitcnt_t shift, mp_bitcnt_t precision);
	void Set_result(mpz_ptr result, mp_size_t size);

	std::vector<uint32_t> x_transform, y_transform, scratch;
	std::vector<uint64_t> coefficients;
	std::vector<mp_limb_t> product, high;
};
