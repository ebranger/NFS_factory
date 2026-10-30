//Integer products with number-theoretic transforms, see int_ntt.h.

#include "int_ntt.h"
#include "arith32.h"
#include <stdlib.h>

//The primes, largest first: 63 2^25 + 1, 15 2^27 + 1, 27 2^26 + 1, 51 2^25 + 1 and 33 2^25 + 1. They are the largest primes
//below 2^31 with 2^25 dividing p - 1, which keeps transforms of 2^25 points possible. A product uses the first P of them.
static const uint32_t ntt_primes[NTT_MAX_PRIMES] = { 2113929217u, 2013265921u, 1811939329u, 1711276033u, 1107296257u };

//The cost of a plan is taken as primes * T * (logT + NTT_COST_OFFSET): the radix-2 levels of the transforms, and the work per
//point that does not depend on the length (splitting, pointwise products, CRT and carries), in units of one level.
#define NTT_COST_OFFSET 6

static uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t p)
{
	uint64_t r = 1;
	a %= p;
	while (e != 0)
	{
		if (e & 1)
		{
			r = r * a % p;
		}
		a = a * a % p;
		e >>= 1;
	}
	return r;
}

static Ntt_tables* build_tables()
{
	Ntt_tables* tables = new Ntt_tables;
	for (int q = 0; q < NTT_MAX_PRIMES; q++)
	{
		const uint32_t p = ntt_primes[q];
		tables->prime[q].Init(p);
		for (int j = 0; j < q; j++)
		{
			tables->garner[j][q] = (uint32_t)pow_mod(ntt_primes[j] % p, p - 2, p);
			tables->garner_shoup[j][q] = shoup_factor(tables->garner[j][q], p);
		}
		for (int logT = 0; logT <= NTT_MAX_LOG_LENGTH; logT++)
		{
			uint64_t s = pow_mod(pow_mod(2, logT, p), p - 2, p) * tables->prime[q].two32 % p;
			tables->scale[q][logT] = (uint32_t)s;
			tables->scale_shoup[q][logT] = shoup_factor((uint32_t)s, p);
		}
	}

	//The bit size of the product of the first P primes, from its 32-bit words.
	uint32_t words[NTT_MAX_PRIMES + 1] = { 1 };
	tables->modulus_bits[0] = 0;
	for (int q = 0; q < NTT_MAX_PRIMES; q++)
	{
		uint64_t carry = 0;
		for (int k = 0; k <= NTT_MAX_PRIMES; k++)
		{
			uint64_t t = (uint64_t)words[k] * ntt_primes[q] + carry;
			words[k] = (uint32_t)t;
			carry = t >> 32;
		}
		int top = NTT_MAX_PRIMES;
		while (words[top] == 0)
		{
			top--;
		}
		int bits = 0;
		while ((words[top] >> bits) > 1)
		{
			bits++;
		}
		tables->modulus_bits[q + 1] = 32 * top + bits; //floor(log2), as the product is not a power of 2
	}
	return tables;
}

const Ntt_tables& ntt_tables()
{
	static const Ntt_tables* tables = build_tables();
	return *tables;
}

static int ceil_log2(size_t n)
{
	int k = 0;
	while (((size_t)1 << k) < n)
	{
		k++;
	}
	return k;
}

//Each coefficient of the convolution is at most n (2^b - 1)^2 < 2^(2b + ceil(log2(n))), which must not exceed the product of the primes.
int ntt_max_coefficient_bits(int primes, size_t n)
{
	int b = (ntt_tables().modulus_bits[primes] - ceil_log2(n)) / 2;
	return (b < NTT_MAX_COEFFICIENT_BITS) ? b : NTT_MAX_COEFFICIENT_BITS;
}

bool ntt_plan_valid(const Ntt_plan& plan, mp_bitcnt_t a_bits, mp_bitcnt_t b_bits, bool cyclic)
{
	const mp_bitcnt_t b = plan.coefficient_bits;
	const size_t T = plan.Length();
	if (plan.primes < 1 || plan.primes > NTT_MAX_PRIMES || plan.log_length < NTT_MIN_LOG_LENGTH || plan.log_length > NTT_MAX_LOG_LENGTH ||
		b < 1 || b > NTT_MAX_COEFFICIENT_BITS)
	{
		return false;
	}
	const mp_bitcnt_t a_count = (a_bits + b - 1) / b, b_count = (b_bits + b - 1) / b;
	if (a_count > T || b_count > T)
	{
		return false;
	}
	if (a_count == 0 || b_count == 0)
	{
		return true;
	}
	if (!cyclic && a_count + b_count - 1 > T)
	{
		return false;
	}
	return (int)b <= ntt_max_coefficient_bits(plan.primes, (size_t)((a_count < b_count) ? a_count : b_count));
}

static double plan_cost(int primes, int logT)
{
	return (double)primes * (double)((size_t)1 << logT) * (logT + NTT_COST_OFFSET);
}

bool ntt_plan_product(Ntt_plan& plan, mp_bitcnt_t a_bits, mp_bitcnt_t b_bits)
{
	double best = -1;
	for (int primes = NTT_MIN_PRIMES; primes <= NTT_MAX_PRIMES; primes++)
	{
		for (int logT = NTT_MIN_LOG_LENGTH; logT <= NTT_MAX_LOG_LENGTH; logT++)
		{
			const mp_bitcnt_t T = (mp_bitcnt_t)1 << logT;
			mp_bitcnt_t b = (a_bits + b_bits + T - 1) / T;
			if (b == 0)
			{
				b = 1;
			}
			while (b <= NTT_MAX_COEFFICIENT_BITS && (a_bits + b - 1) / b + (b_bits + b - 1) / b - 1 > T)
			{
				b++;
			}
			Ntt_plan candidate = { primes, logT, (int)b };
			if (b <= NTT_MAX_COEFFICIENT_BITS && ntt_plan_valid(candidate, a_bits, b_bits, false) && (best < 0 || plan_cost(primes, logT) < best))
			{
				best = plan_cost(primes, logT);
				plan = candidate;
			}
		}
	}
	return best >= 0;
}

bool ntt_plan_cyclic(Ntt_plan& plan, mp_bitcnt_t min_L, mp_bitcnt_t a_bits, mp_bitcnt_t b_bits)
{
	mp_bitcnt_t L = min_L;
	L = (a_bits > L) ? a_bits : L;
	L = (b_bits > L) ? b_bits : L;
	double best = -1;
	for (int primes = NTT_MIN_PRIMES; primes <= NTT_MAX_PRIMES; primes++)
	{
		for (int logT = NTT_MIN_LOG_LENGTH; logT <= NTT_MAX_LOG_LENGTH; logT++)
		{
			const mp_bitcnt_t T = (mp_bitcnt_t)1 << logT;
			mp_bitcnt_t b = (L + T - 1) / T;
			if (b == 0)
			{
				b = 1;
			}
			Ntt_plan candidate = { primes, logT, (int)((b <= NTT_MAX_COEFFICIENT_BITS) ? b : 0) };
			if (ntt_plan_valid(candidate, a_bits, b_bits, true) && (best < 0 || plan_cost(primes, logT) < best))
			{
				best = plan_cost(primes, logT);
				plan = candidate;
			}
		}
	}
	return best >= 0;
}

//The bits of a from bit pos upwards, b of them (b <= 64). Bit pos must be below the top of a.
static inline uint64_t get_bits(mp_srcptr a, mp_size_t a_limbs, mp_bitcnt_t pos, int b)
{
	const mp_size_t i = (mp_size_t)(pos / GMP_LIMB_BITS);
	const unsigned int shift = (unsigned int)(pos % GMP_LIMB_BITS);
	uint64_t v = a[i] >> shift;
	if (shift != 0 && i + 1 < a_limbs)
	{
		v |= a[i + 1] << (GMP_LIMB_BITS - shift);
	}
	return (b < 64) ? v & ((1ULL << b) - 1) : v;
}

//x mod p for a 32-bit x, with p > 2^30.
static inline uint32_t reduce_word(uint32_t x, uint32_t p)
{
	while (x >= p)
	{
		x -= p;
	}
	return x;
}

void Ntt_product::Forward(std::vector<uint32_t>& transform, const Ntt_plan& plan, mp_srcptr a, mp_size_t a_limbs)
{
	const Ntt_tables& tables = ntt_tables();
	const size_t T = plan.Length();
	const int b = plan.coefficient_bits;

	while (a_limbs > 0 && a[a_limbs - 1] == 0)
	{
		a_limbs--;
	}
	size_t count = 0;
	if (a_limbs > 0)
	{
		mp_bitcnt_t bits = (mp_bitcnt_t)a_limbs * GMP_LIMB_BITS;
		for (mp_limb_t top = a[a_limbs - 1]; (top >> (GMP_LIMB_BITS - 1)) == 0; top <<= 1)
		{
			bits--;
		}
		count = (size_t)((bits + b - 1) / b);
	}
	coefficients.resize(count);
	for (size_t i = 0; i < count; i++)
	{
		coefficients[i] = get_bits(a, a_limbs, (mp_bitcnt_t)i * b, b);
	}

	transform.resize(plan.primes * T);
	if (scratch.size() < ntt_scratch_size(plan.log_length))
	{
		scratch.resize(ntt_scratch_size(plan.log_length));
	}
	for (int q = 0; q < plan.primes; q++)
	{
		const Ntt_prime& prime = tables.prime[q];
		const uint32_t p = prime.p;
		uint32_t* x = transform.data() + q * T;
		for (size_t i = 0; i < count; i++)
		{
			uint64_t v = coefficients[i];
			x[i] = add_mod(mul_shoup((uint32_t)(v >> 32), prime.two32, prime.two32_shoup, p), reduce_word((uint32_t)v, p), p);
		}
		for (size_t i = count; i < T; i++)
		{
			x[i] = 0;
		}
		ntt_forward(x, plan.log_length, prime, scratch.data());
	}
}

void Ntt_product::Pointwise(std::vector<uint32_t>& x, const std::vector<uint32_t>& y, const Ntt_plan& plan)
{
	const Ntt_tables& tables = ntt_tables();
	const size_t T = plan.Length();
	for (int q = 0; q < plan.primes; q++)
	{
		ntt_pointwise(x.data() + q * T, y.data() + q * T, T, tables.prime[q]);
	}
}

//*x += y + carry, returning the carry out.
static inline unsigned int add_carry(mp_limb_t* x, mp_limb_t y, unsigned int carry)
{
	mp_limb_t s = *x + y;
	unsigned int c = (s < y);
	mp_limb_t t = s + carry;
	c |= (t < s);
	*x = t;
	return c;
}

mp_size_t Ntt_product::Inverse(std::vector<mp_limb_t>& result, std::vector<uint32_t>& transform, const Ntt_plan& plan, bool cyclic)
{
	const Ntt_tables& tables = ntt_tables();
	const int P = plan.primes, logT = plan.log_length;
	const size_t T = plan.Length();
	const mp_bitcnt_t b = plan.coefficient_bits;

	if (scratch.size() < ntt_scratch_size(logT))
	{
		scratch.resize(ntt_scratch_size(logT));
	}
	for (int q = 0; q < P; q++)
	{
		ntt_inverse(transform.data() + q * T, logT, tables.prime[q], scratch.data());
	}

	//Each coefficient is less than the product of the primes, below 2^160, and is added at bit k b. The sum is less than
	//2^((T + 1) b + logT + 1).
	const mp_size_t limbs = (mp_size_t)((T - 1) * b / GMP_LIMB_BITS) + 5;
	result.assign(limbs, 0);
	mp_ptr r = result.data();
	uint32_t p[NTT_MAX_PRIMES], scale[NTT_MAX_PRIMES], scale_shoup[NTT_MAX_PRIMES];
	for (int q = 0; q < P; q++)
	{
		p[q] = tables.prime[q].p;
		scale[q] = tables.scale[q][logT];
		scale_shoup[q] = tables.scale_shoup[q][logT];
	}

	for (size_t k = 0; k < T; k++)
	{
		//The residues of the coefficient, without the factor T of the inverse transform and the 1 / 2^32 of the pointwise product.
		uint32_t v[NTT_MAX_PRIMES];
		for (int q = 0; q < P; q++)
		{
			v[q] = mul_shoup(transform[q * T + k], scale[q], scale_shoup[q], p[q]);
		}

		//Garner's mixed radix CRT: c = v0 + p0 (v1 + p1 (v2 + ...)).
		for (int i = 1; i < P; i++)
		{
			uint32_t t = v[i];
			for (int j = 0; j < i; j++)
			{
				uint32_t vj = (v[j] >= p[i]) ? v[j] - p[i] : v[j]; //v[j] < p[j] < 2^31 < 2 p[i]
				t = mul_shoup(sub_mod(t, vj, p[i]), tables.garner[j][i], tables.garner_shoup[j][i], p[i]);
			}
			v[i] = t;
		}
		uint32_t words[6] = { v[P - 1], 0, 0, 0, 0, 0 };
		for (int i = P - 2; i >= 0; i--)
		{
			uint64_t carry = v[i];
			for (int w = 0; w < P; w++)
			{
				uint64_t t = (uint64_t)words[w] * p[i] + carry;
				words[w] = (uint32_t)t;
				carry = t >> 32;
			}
		}
		mp_limb_t c[3] = { words[0] | ((mp_limb_t)words[1] << 32), words[2] | ((mp_limb_t)words[3] << 32), words[4] };

		//Add the coefficient at bit k b.
		const mp_bitcnt_t pos = (mp_bitcnt_t)k * b;
		mp_ptr dst = r + pos / GMP_LIMB_BITS;
		const unsigned int shift = (unsigned int)(pos % GMP_LIMB_BITS);
		mp_limb_t s[4];
		if (shift == 0)
		{
			s[0] = c[0];
			s[1] = c[1];
			s[2] = c[2];
			s[3] = 0;
		}
		else
		{
			s[0] = c[0] << shift;
			s[1] = (c[1] << shift) | (c[0] >> (GMP_LIMB_BITS - shift));
			s[2] = (c[2] << shift) | (c[1] >> (GMP_LIMB_BITS - shift));
			s[3] = c[2] >> (GMP_LIMB_BITS - shift);
		}
		unsigned int carry = 0;
		for (int w = 0; w < 4; w++)
		{
			carry = add_carry(dst + w, s[w], carry);
		}
		for (mp_ptr d = dst + 4; carry != 0; d++)
		{
			carry = add_carry(d, 0, carry);
		}
	}

	if (!cyclic)
	{
		mp_size_t size = limbs;
		while (size > 0 && r[size - 1] == 0)
		{
			size--;
		}
		return size;
	}

	//Modulo 2^L - 1: add the part from bit L upwards onto the bottom. It is less than 2^L, so the sum F is less than 2^(L + 1).
	const mp_bitcnt_t L = plan.L();
	const mp_size_t L_limbs = (mp_size_t)(L / GMP_LIMB_BITS);
	const unsigned int L_bit = (unsigned int)(L % GMP_LIMB_BITS);
	const mp_limb_t L_mask = (((mp_limb_t)1) << L_bit) - 1;
	mp_size_t high_limbs = limbs - L_limbs;
	high.resize(high_limbs);
	if (L_bit != 0)
	{
		mpn_rshift(high.data(), r + L_limbs, high_limbs, L_bit);
		r[L_limbs] &= L_mask;
		mpn_zero(r + L_limbs + 1, limbs - L_limbs - 1);
	}
	else
	{
		mpn_copyi(high.data(), r + L_limbs, high_limbs);
		mpn_zero(r + L_limbs, high_limbs);
	}
	while (high_limbs > 0 && high[high_limbs - 1] == 0)
	{
		high_limbs--;
	}
	const mp_size_t size = L_limbs + 1;
	if (high_limbs > 0)
	{
		mpn_add(r, r, size, high.data(), high_limbs);
	}

	//The end-around carry: 2^L is 1. Then F is at most 2^L - 1, which is 0.
	if ((r[L_limbs] >> L_bit) & 1)
	{
		r[L_limbs] &= L_mask;
		mpn_add_1(r, r, size, 1);
	}
	bool all_ones = (L_bit == 0 || r[L_limbs] == L_mask);
	for (mp_size_t i = 0; all_ones && i < L_limbs; i++)
	{
		all_ones = (r[i] == ~(mp_limb_t)0);
	}
	if (all_ones)
	{
		mpn_zero(r, size);
	}

	mp_size_t n = size;
	while (n > 0 && r[n - 1] == 0)
	{
		n--;
	}
	return n;
}

void Ntt_product::Set_result(mpz_ptr result, mp_size_t size)
{
	if (result->_mp_alloc < size)
	{
		mpz_realloc2(result, size * GMP_LIMB_BITS);
	}
	mpn_copyi(result->_mp_d, product.data(), size);
	result->_mp_size = (int)size;
}

void Ntt_product::Forward(std::vector<uint32_t>& transform, const Ntt_plan& plan, mpz_srcptr a)
{
	Forward(transform, plan, a->_mp_d, abs(a->_mp_size));
}

void Ntt_product::Product_from(mpz_ptr result, std::vector<uint32_t>& transform, const Ntt_plan& plan)
{
	Set_result(result, Inverse(product, transform, plan, false));
}

void Ntt_product::Middle_from(mpz_ptr result, std::vector<uint32_t>& transform, const Ntt_plan& plan, mp_bitcnt_t shift, mp_bitcnt_t precision)
{
	Inverse(product, transform, plan, true);
	Extract_middle(result, shift, precision);
}

void Ntt_product::Multiply(mpz_ptr result, mpz_srcptr a, mpz_srcptr b)
{
	const mp_size_t a_limbs = abs(a->_mp_size), b_limbs = abs(b->_mp_size);
	if (a_limbs == 0 || b_limbs == 0)
	{
		mpz_set_ui(result, 0);
		return;
	}
	Ntt_plan plan;
	if (!ntt_plan_product(plan, mpz_sizeinbase(a, 2), mpz_sizeinbase(b, 2)))
	{
		mpz_mul(result, a, b); //Beyond the longest transform.
		return;
	}
	const bool negative = (a->_mp_size < 0) != (b->_mp_size < 0);

	Forward(x_transform, plan, a->_mp_d, a_limbs);
	if (a == b)
	{
		Pointwise(x_transform, x_transform, plan);
	}
	else
	{
		Forward(y_transform, plan, b->_mp_d, b_limbs);
		Pointwise(x_transform, y_transform, plan);
	}
	Set_result(result, Inverse(product, x_transform, plan, false));
	if (negative)
	{
		mpz_neg(result, result);
	}
}

mp_bitcnt_t Ntt_product::Product_mod(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t min_L)
{
	const mp_size_t y_limbs = abs(y->_mp_size), x_limbs = abs(x->_mp_size);
	const mp_bitcnt_t y_bits = mpz_sizeinbase(y, 2), x_bits = mpz_sizeinbase(x, 2);
	Ntt_plan plan;
	if (!ntt_plan_cyclic(plan, min_L, y_bits, x_bits))
	{
		//Beyond the longest transform: the exact product, reduced modulo 2^L - 1 with L = min_L.
		mpz_t high_part;
		mpz_init(high_part);
		mpz_mul(result, y, x);
		while (mpz_sizeinbase(result, 2) > min_L)
		{
			mpz_tdiv_q_2exp(high_part, result, min_L);
			mpz_tdiv_r_2exp(result, result, min_L);
			mpz_add(result, result, high_part);
		}
		mpz_add_ui(result, result, 1);
		if (mpz_sizeinbase(result, 2) > min_L)
		{
			mpz_set_ui(result, 0);
		}
		else
		{
			mpz_sub_ui(result, result, 1);
		}
		mpz_clear(high_part);
		return min_L;
	}
	if (y_limbs == 0 || x_limbs == 0)
	{
		mpz_set_ui(result, 0);
		return plan.L();
	}

	Forward(x_transform, plan, x->_mp_d, x_limbs);
	Forward(y_transform, plan, y->_mp_d, y_limbs);
	Pointwise(x_transform, y_transform, plan);
	Set_result(result, Inverse(product, x_transform, plan, true));
	return plan.L();
}

//From the product modulo 2^L - 1, result = bits from shift up to precision, minus one, as in Cyclic_product.
void Ntt_product::Extract_middle(mpz_ptr result, mp_bitcnt_t shift, mp_bitcnt_t precision)
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

//The product modulo 2^L - 1 needs L at least 64 bits more than precision + bits of x - shift, as in Cyclic_product, and at least
//64 bits more than precision. The part of the product above 2^L is then less than 2^(shift - 64). It wraps around and is added to the
//low bits, below the bits that are needed, which can make the result one too large; so can taking 2^L - 1 as 0. One is subtracted,
//so that the result is correct or one too small.
static mp_bitcnt_t middle_min_L(mp_bitcnt_t x_bits, mp_bitcnt_t shift, mp_bitcnt_t precision)
{
	return ((x_bits > shift) ? precision + x_bits - shift : precision) + GMP_LIMB_BITS;
}

void Ntt_product::Middle_product(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t shift, mp_bitcnt_t precision)
{
	const mp_size_t y_limbs = abs(y->_mp_size), x_limbs = abs(x->_mp_size);
	if (y_limbs == 0 || x_limbs == 0)
	{
		mpz_set_ui(result, 0);
		return;
	}

	const mp_bitcnt_t x_bits = mpz_sizeinbase(x, 2);
	Ntt_plan plan;
	if (!ntt_plan_cyclic(plan, middle_min_L(x_bits, shift, precision), mpz_sizeinbase(y, 2), x_bits))
	{
		mpz_mul(result, y, x); //Beyond the longest transform: exact.
		mpz_tdiv_r_2exp(result, result, precision);
		mpz_tdiv_q_2exp(result, result, shift);
		return;
	}
	Forward(y_transform, plan, y->_mp_d, y_limbs);
	Forward(x_transform, plan, x->_mp_d, x_limbs);
	Pointwise(x_transform, y_transform, plan);
	Inverse(product, x_transform, plan, true);
	Extract_middle(result, shift, precision);
}

void Ntt_product::Middle_product_pair(mpz_ptr result0, mpz_ptr result1, mpz_srcptr y, mpz_srcptr x0, mpz_srcptr x1, mp_bitcnt_t shift0, mp_bitcnt_t shift1, mp_bitcnt_t precision)
{
	const mp_size_t y_limbs = abs(y->_mp_size), x0_limbs = abs(x0->_mp_size), x1_limbs = abs(x1->_mp_size);
	if (y_limbs == 0 || x0_limbs == 0 || x1_limbs == 0)
	{
		Middle_product(result0, y, x0, shift0, precision);
		Middle_product(result1, y, x1, shift1, precision);
		return;
	}

	const mp_bitcnt_t x0_bits = mpz_sizeinbase(x0, 2), x1_bits = mpz_sizeinbase(x1, 2);
	const mp_bitcnt_t min_L0 = middle_min_L(x0_bits, shift0, precision), min_L1 = middle_min_L(x1_bits, shift1, precision);
	Ntt_plan plan;
	if (!ntt_plan_cyclic(plan, (min_L0 > min_L1) ? min_L0 : min_L1, mpz_sizeinbase(y, 2), (x0_bits > x1_bits) ? x0_bits : x1_bits))
	{
		Middle_product(result0, y, x0, shift0, precision);
		Middle_product(result1, y, x1, shift1, precision);
		return;
	}

	Forward(y_transform, plan, y->_mp_d, y_limbs);

	Forward(x_transform, plan, x0->_mp_d, x0_limbs);
	Pointwise(x_transform, y_transform, plan);
	Inverse(product, x_transform, plan, true);
	Extract_middle(result0, shift0, precision);

	Forward(x_transform, plan, x1->_mp_d, x1_limbs);
	Pointwise(x_transform, y_transform, plan);
	Inverse(product, x_transform, plan, true);
	Extract_middle(result1, shift1, precision);
}
