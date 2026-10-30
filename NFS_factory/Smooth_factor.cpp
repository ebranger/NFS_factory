//Factoring of smooth relation values with native 64-bit arithmetic.
//After the large primes are divided out, the value of a relation from the batch check only has prime factors below
//the factor base limit, which is less than 2^32. The value then fits in 128 bits, and any remaining cofactor of at
//least 2^32 must be composite. Small primes are found by trial division, the rest by Pollard-Brent rho with Montgomery
//multiplication modulo the (up to 128-bit) cofactor.
//The same arithmetic gives a fast primality test for the cofactors of the batch check.

#include "Smooth_factor.h"
#include <intrin.h>
#include <algorithm>
#include <vector>

typedef unsigned long long u64;
typedef unsigned int u32;

#define RHO_MAX_ITERATIONS (1 << 20) //sqrt(2^32) = 65536 iterations are expected for the largest primes, so this is plenty.
#define RHO_MPZ_MAX_ITERATIONS (1 << 19) //The same for rho in MPIR arithmetic on values of more than 128 bits, which is slower, so it gives up sooner.
#define RHO_GCD_INTERVAL 64 //Iterations between the gcd computations in rho.
#define MAX_PENDING 128 //Most cofactors waiting to be split.

static const u32 small_primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };

//A 128-bit unsigned number as two 64-bit limbs.
struct u128
{
	u64 lo, hi;
};

static inline bool is_zero(const u128& a)
{
	return (a.lo | a.hi) == 0;
}

static inline bool equal(const u128& a, const u128& b)
{
	return a.lo == b.lo && a.hi == b.hi;
}

static inline bool less(const u128& a, const u128& b)
{
	return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}

static inline u128 subtract(const u128& a, const u128& b)
{
	u128 r;
	unsigned char borrow = _subborrow_u64(0, a.lo, b.lo, &r.lo);
	_subborrow_u64(borrow, a.hi, b.hi, &r.hi);
	return r;
}

static inline u128 shift_right(const u128& a, unsigned int bits)
{
	u128 r;
	if (bits == 0)
	{
		return a;
	}
	if (bits >= 64)
	{
		r.lo = a.hi >> (bits - 64);
		r.hi = 0;
	}
	else
	{
		r.lo = (a.lo >> bits) | (a.hi << (64 - bits));
		r.hi = a.hi >> bits;
	}
	return r;
}

//Divide a by d, with d < 2^64. Returns the quotient and stores the remainder.
static inline u128 divide(const u128& a, u64 d, u64* remainder)
{
	u128 q;
	q.hi = a.hi / d;
	q.lo = _udiv128(a.hi % d, a.lo, d, remainder);
	return q;
}

//Greatest common divisor of a and an odd b, using the binary algorithm.
static u128 gcd(u128 a, u128 b)
{
	while (!is_zero(a))
	{
		unsigned long zeros;
		if (a.lo != 0)
		{
			_BitScanForward64(&zeros, a.lo);
		}
		else
		{
			_BitScanForward64(&zeros, a.hi);
			zeros += 64;
		}
		a = shift_right(a, zeros);

		//Both are odd now.
		if (less(a, b))
		{
			u128 t = a;
			a = b;
			b = t;
		}
		a = subtract(a, b);
	}
	return b;
}

//-1/n mod 2^64 for odd n, by Newton iteration.
static u64 montgomery_inverse(u64 n)
{
	u64 inverse = n; //correct to 3 bits, since n * n = 1 mod 8 for odd n.
	for (int i = 0; i < 5; i++)
	{
		inverse *= 2 - n * inverse;
	}
	return 0 - inverse;
}

//r + carry * 2^128 mod n, for r + carry * 2^128 < 2n. Without branches, since whether n is subtracted is unpredictable.
static inline u128 reduce_once(const u128& r, u64 carry, const u128& n)
{
	u128 s;
	unsigned char borrow = _subborrow_u64(0, r.lo, n.lo, &s.lo);
	borrow = _subborrow_u64(borrow, r.hi, n.hi, &s.hi);
	u64 keep = 0 - (u64)(borrow & (carry == 0)); //all ones if r < n.
	s.lo ^= (s.lo ^ r.lo) & keep;
	s.hi ^= (s.hi ^ r.hi) & keep;
	return s;
}

//Montgomery product a * b / 2^128 mod n, for a, b < n and odd n. The result is less than n.
static inline u128 montgomery_multiply(const u128& a, const u128& b, const u128& n, u64 ninv)
{
	u64 t0 = 0, t1 = 0, t2 = 0, t3;
	u64 b_limbs[2] = { b.lo, b.hi };

	for (int i = 0; i < 2; i++)
	{
		u64 hi, lo, carry;
		unsigned char c;

		//t += a * b_i
		lo = _umul128(a.lo, b_limbs[i], &hi);
		c = _addcarry_u64(0, t0, lo, &t0);
		carry = hi + c;
		lo = _umul128(a.hi, b_limbs[i], &hi);
		c = _addcarry_u64(0, t1, lo, &t1);
		hi += c;
		c = _addcarry_u64(0, t1, carry, &t1);
		hi += c;
		c = _addcarry_u64(0, t2, hi, &t2);
		t3 = c;

		//t = (t + m * n) / 2^64, where m makes the lowest limb zero.
		u64 m = t0 * ninv;
		lo = _umul128(m, n.lo, &hi);
		c = _addcarry_u64(0, t0, lo, &t0);
		carry = hi + c;
		lo = _umul128(m, n.hi, &hi);
		c = _addcarry_u64(0, t1, lo, &t1);
		hi += c;
		c = _addcarry_u64(0, t1, carry, &t1);
		hi += c;
		c = _addcarry_u64(0, t2, hi, &t2);
		t3 += c;

		t0 = t1;
		t1 = t2;
		t2 = t3;
	}

	u128 r = { t0, t1 };
	return reduce_once(r, t2, n);
}

//Montgomery square a^2 / 2^128 mod n, for a < n and odd n. The same as montgomery_multiply(a, a, n, ninv), unrolled,
//and with the cross product computed once.
static inline u128 montgomery_square(const u128& a, const u128& n, u64 ninv)
{
	u64 t0, t1, t2, t3, hi, lo, h, c4 = 0;
	unsigned char c;

	//t = a^2 = a.lo^2 + 2 a.lo a.hi 2^64 + a.hi^2 2^128. It is less than 2^256, so the last carry is 0.
	t0 = _umul128(a.lo, a.lo, &t1);
	t2 = _umul128(a.hi, a.hi, &t3);
	lo = _umul128(a.lo, a.hi, &hi);
	u64 top = hi >> 63;
	hi = (hi << 1) | (lo >> 63);
	lo <<= 1;
	c = _addcarry_u64(0, t1, lo, &t1);
	c = _addcarry_u64(c, t2, hi, &t2);
	_addcarry_u64(c, t3, top, &t3);

	//Two Montgomery reduction steps, each clearing the lowest limb. c4 is the carry out of t3.
	u64 m = t0 * ninv;
	lo = _umul128(m, n.lo, &hi);
	c = _addcarry_u64(0, t0, lo, &t0);
	h = hi;
	lo = _umul128(m, n.hi, &hi);
	c = _addcarry_u64(c, t1, h, &t1);
	hi += c;
	c = _addcarry_u64(0, t1, lo, &t1);
	c = _addcarry_u64(c, t2, hi, &t2);
	c = _addcarry_u64(c, t3, 0, &t3);
	c4 += c;

	m = t1 * ninv;
	lo = _umul128(m, n.lo, &hi);
	c = _addcarry_u64(0, t1, lo, &t1);
	h = hi;
	lo = _umul128(m, n.hi, &hi);
	c = _addcarry_u64(c, t2, h, &t2);
	hi += c;
	c = _addcarry_u64(0, t2, lo, &t2);
	c = _addcarry_u64(c, t3, hi, &t3);
	c4 += c;

	u128 r = { t2, t3 };
	return reduce_once(r, c4, n);
}

//The rho map y -> y^2 + c mod n. Working on Montgomery representations is fine, since the result mod any prime p
//dividing n only depends on y mod p.
static inline u128 rho_step(const u128& y, const u128& c, const u128& n, u64 ninv)
{
	u128 r = montgomery_square(y, n, ninv);
	u128 sum;
	unsigned char carry = _addcarry_u64(0, r.lo, c.lo, &sum.lo);
	carry = _addcarry_u64(carry, r.hi, c.hi, &sum.hi);
	return reduce_once(sum, carry, n);
}

//|a - b|, without branches.
static inline u128 difference(const u128& a, const u128& b)
{
	u128 d, e;
	unsigned char borrow = _subborrow_u64(0, a.lo, b.lo, &d.lo);
	borrow = _subborrow_u64(borrow, a.hi, b.hi, &d.hi);
	u64 negative = 0 - (u64)borrow;
	_subborrow_u64(_subborrow_u64(0, 0, d.lo, &e.lo), 0, d.hi, &e.hi); //e = -d
	d.lo ^= (d.lo ^ e.lo) & negative;
	d.hi ^= (d.hi ^ e.hi) & negative;
	return d;
}

//Find a nontrivial factor of the odd composite n with Pollard-Brent rho. Returns false if none was found.
static bool rho(const u128& n, u128* factor)
{
	u64 ninv = montgomery_inverse(n.lo);
	u128 one = { 1, 0 };

	for (u64 constant = 1; constant <= 5; constant++)
	{
		u128 c = { constant, 0 };
		u128 y = { 2, 0 }, x = y, saved_y = y, q = one, g = one;
		u64 r = 1;

		while (equal(g, one) && r <= RHO_MAX_ITERATIONS)
		{
			x = y;
			for (u64 i = 0; i < r; i++)
			{
				y = rho_step(y, c, n, ninv);
			}

			for (u64 k = 0; k < r && equal(g, one); k += RHO_GCD_INTERVAL)
			{
				saved_y = y;
				u64 steps = std::min((u64)RHO_GCD_INTERVAL, r - k);
				for (u64 i = 0; i < steps; i++)
				{
					y = rho_step(y, c, n, ninv);
					q = montgomery_multiply(q, difference(x, y), n, ninv);
				}
				g = gcd(q, n);
			}
			r *= 2;
		}

		if (equal(g, n))
		{
			//Several factors were found at once, go back and take one step at a time.
			do
			{
				saved_y = rho_step(saved_y, c, n, ninv);
				g = gcd(difference(x, saved_y), n);
			} while (equal(g, one));
		}

		if (!equal(g, one) && !equal(g, n))
		{
			*factor = g;
			return true;
		}
	}
	return false;
}

//Montgomery product a * b / 2^64 mod n, for a, b < n < 2^64 and odd n. The result is less than n.
static inline u64 montgomery_multiply64(u64 a, u64 b, u64 n, u64 ninv)
{
	u64 hi, lo = _umul128(a, b, &hi);
	u64 m = lo * ninv;
	u64 mhi, mlo = _umul128(m, n, &mhi);

	//(lo + mlo) is 0 mod 2^64, with a carry unless lo is 0.
	u64 t;
	unsigned char carry = _addcarry_u64(lo != 0, hi, mhi, &t);
	if (carry || t >= n)
	{
		t -= n;
	}
	return t;
}

static inline u64 rho_step64(u64 y, u64 c, u64 n, u64 ninv)
{
	u64 r = montgomery_multiply64(y, y, n, ninv);
	u64 sum;
	unsigned char carry = _addcarry_u64(0, r, c, &sum);
	if (carry || sum >= n)
	{
		sum -= n;
	}
	return sum;
}

static u64 gcd64(u64 a, u64 b)
{
	while (b != 0)
	{
		u64 t = a % b;
		a = b;
		b = t;
	}
	return a;
}

//Pollard-Brent rho for an odd composite n < 2^64. Returns false if no factor was found.
static bool rho64(u64 n, u64* factor)
{
	u64 ninv = montgomery_inverse(n);

	for (u64 c = 1; c <= 5; c++)
	{
		u64 y = 2, x = y, saved_y = y, q = 1, g = 1;
		u64 r = 1;

		while (g == 1 && r <= RHO_MAX_ITERATIONS)
		{
			x = y;
			for (u64 i = 0; i < r; i++)
			{
				y = rho_step64(y, c, n, ninv);
			}

			for (u64 k = 0; k < r && g == 1; k += RHO_GCD_INTERVAL)
			{
				saved_y = y;
				u64 steps = std::min((u64)RHO_GCD_INTERVAL, r - k);
				for (u64 i = 0; i < steps; i++)
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
			*factor = g;
			return true;
		}
	}
	return false;
}

static u64 powmod32(u64 base, u64 exponent, u64 n)
{
	u64 result = 1;
	base %= n;
	while (exponent > 0)
	{
		if (exponent & 1)
		{
			result = result * base % n;
		}
		base = base * base % n;
		exponent >>= 1;
	}
	return result;
}

//Deterministic Miller-Rabin test for n < 2^32, with the bases 2, 7 and 61.
static bool is_prime32(u64 n)
{
	if (n < 2)
	{
		return false;
	}
	static const u64 bases[] = { 2, 7, 61 };
	for (u64 p : bases)
	{
		if (n == p)
		{
			return true;
		}
		if (n % p == 0)
		{
			return false;
		}
	}

	u64 d = n - 1;
	int s = 0;
	while ((d & 1) == 0)
	{
		d >>= 1;
		s++;
	}

	for (u64 a : bases)
	{
		u64 x = powmod32(a, d, n);
		if (x == 1 || x == n - 1)
		{
			continue;
		}
		bool composite = true;
		for (int i = 1; i < s; i++)
		{
			x = x * x % n;
			if (x == n - 1)
			{
				composite = false;
				break;
			}
		}
		if (composite)
		{
			return false;
		}
	}
	return true;
}

int factor_smooth(unsigned int* factors, mpz_srcptr value, int max_factors)
{
	int size = abs(value->_mp_size);
	if (size == 0 || size > 2)
	{
		return -1;
	}

	u128 n = { value->_mp_d[0], (size > 1) ? value->_mp_d[1] : 0 };
	u128 one = { 1, 0 };
	int count = 0;

	//Factors of 2.
	while ((n.lo & 1) == 0)
	{
		if (count >= max_factors)
		{
			return -1;
		}
		factors[count++] = 2;
		n = shift_right(n, 1);
	}

	//Trial division by the smallest odd primes.
	for (u32 p : small_primes)
	{
		while (true)
		{
			u64 remainder;
			u128 q = divide(n, p, &remainder);
			if (remainder != 0)
			{
				break;
			}
			if (count >= max_factors)
			{
				return -1;
			}
			factors[count++] = p;
			n = q;
		}
	}

	//Split the rest. All cofactors are odd, and have no prime factors below 100.
	u128 pending[MAX_PENDING];
	int num_pending = 0;
	if (!equal(n, one))
	{
		pending[num_pending++] = n;
	}

	while (num_pending > 0)
	{
		u128 c = pending[--num_pending];

		if (c.hi == 0 && c.lo < ((u64)1 << 32) && is_prime32(c.lo))
		{
			if (count >= max_factors)
			{
				return -1;
			}
			factors[count++] = (u32)c.lo;
			continue;
		}

		//c is composite, unless it has a prime factor above 2^32, in which case rho fails.
		u128 d = { 0, 0 };
		if (num_pending + 2 > MAX_PENDING)
		{
			return -1;
		}

		if (c.hi == 0)
		{
			//Rho with 64-bit arithmetic was about twice as fast as SQUFOF here.
			if (!rho64(c.lo, &d.lo))
			{
				return -1;
			}
		}
		else if (!rho(c, &d))
		{
			return -1;
		}

		u128 e = { 0, 0 };
		if (d.hi != 0)
		{
			//d is a product of several primes, of 2^64 or more, so the cofactor e = c / d is less than 2^64. Since d is odd and
			//c = d e, e is c times the inverse of d modulo 2^64.
			e.lo = c.lo * (0 - montgomery_inverse(d.lo));
		}
		else
		{
			u64 remainder;
			e = divide(c, d.lo, &remainder);
		}
		pending[num_pending++] = d;
		pending[num_pending++] = e;
	}

	std::sort(factors, factors + count);
	return count;
}

//Find a nontrivial factor of the odd composite n, of more than 128 bits, with Pollard-Brent rho in MPIR arithmetic.
//Returns false if none was found.
static bool rho_mpz(mpz_srcptr n, mpz_ptr factor)
{
	mpz_t x, y, saved_y, q, t;
	mpz_inits(x, y, saved_y, q, t, NULL);
	bool found = false;

	for (unsigned long c = 1; c <= 3 && !found; c++)
	{
		mpz_set_ui(y, 2);
		mpz_set_ui(q, 1);
		mpz_set_ui(factor, 1);
		u64 r = 1;

		while (mpz_cmp_ui(factor, 1) == 0 && r <= RHO_MPZ_MAX_ITERATIONS)
		{
			mpz_set(x, y);
			for (u64 i = 0; i < r; i++)
			{
				mpz_mul(y, y, y);
				mpz_add_ui(y, y, c);
				mpz_mod(y, y, n);
			}

			for (u64 k = 0; k < r && mpz_cmp_ui(factor, 1) == 0; k += RHO_GCD_INTERVAL)
			{
				mpz_set(saved_y, y);
				u64 steps = std::min((u64)RHO_GCD_INTERVAL, r - k);
				for (u64 i = 0; i < steps; i++)
				{
					mpz_mul(y, y, y);
					mpz_add_ui(y, y, c);
					mpz_mod(y, y, n);
					mpz_sub(t, x, y);
					mpz_mul(q, q, t);
					mpz_mod(q, q, n);
				}
				mpz_gcd(factor, q, n);
			}
			r *= 2;
		}

		if (mpz_cmp(factor, n) == 0 || mpz_sgn(factor) == 0)
		{
			//Several factors were found at once, go back and take one step at a time.
			do
			{
				mpz_mul(saved_y, saved_y, saved_y);
				mpz_add_ui(saved_y, saved_y, c);
				mpz_mod(saved_y, saved_y, n);
				mpz_sub(t, x, saved_y);
				mpz_gcd(factor, t, n);
			} while (mpz_cmp_ui(factor, 1) == 0);
		}

		found = mpz_cmp_ui(factor, 1) != 0 && mpz_cmp(factor, n) != 0;
	}

	mpz_clears(x, y, saved_y, q, t, NULL);
	return found;
}

//The odd primes below 2^16, for trial division of values of more than 128 bits.
static const std::vector<u32>& trial_division_primes()
{
	static const std::vector<u32> primes = []
	{
		const u32 limit = 1 << 16;
		std::vector<bool> composite(limit, false);
		std::vector<u32> list;
		for (u32 i = 3; i < limit; i += 2)
		{
			if (composite[i])
			{
				continue;
			}
			list.push_back(i);
			for (u32 j = i * i; j < limit; j += 2 * i)
			{
				composite[j] = true;
			}
		}
		return list;
	}();
	return primes;
}

int factor_general(unsigned int* factors, mpz_srcptr value, int max_factors)
{
	if (abs(value->_mp_size) <= 2)
	{
		int count = factor_smooth(factors, value, max_factors);
		if (count >= 0 || value->_mp_size == 0)
		{
			return count;
		}
		//Rho can fail on values with many small primes; trial division below takes most of them out first.
	}

	mpz_t n;
	mpz_init(n);
	mpz_abs(n, value);
	int count = 0;
	bool failed = false;

	//Factors of 2, and trial division by the odd primes below 2^16. Rho alone does badly on values with many small primes: their
	//cycles are short, and they tend to be found all at once, which gives no factor.
	unsigned long twos = (unsigned long)mpz_scan1(n, 0);
	for (unsigned long i = 0; i < twos && !failed; i++)
	{
		failed = count >= max_factors;
		if (!failed)
		{
			factors[count++] = 2;
		}
	}
	mpz_tdiv_q_2exp(n, n, twos);
	for (u32 p : trial_division_primes())
	{
		while (!failed && mpz_divisible_ui_p(n, p))
		{
			failed = count >= max_factors;
			if (!failed)
			{
				factors[count++] = p;
				mpz_divexact_ui(n, n, p);
			}
		}
	}

	//Split the rest until each part has at most 128 bits, and factor those natively.
	mpz_t pending[MAX_PENDING], d;
	int num_pending = 0;
	mpz_init(d);
	if (mpz_cmp_ui(n, 1) != 0)
	{
		mpz_init_set(pending[num_pending++], n);
	}

	while (num_pending > 0 && !failed)
	{
		mpz_ptr c = pending[num_pending - 1];

		if (abs(c->_mp_size) <= 2)
		{
			int found = factor_smooth(factors + count, c, max_factors - count);
			failed = found < 0;
			count += failed ? 0 : found;
			mpz_clear(pending[--num_pending]);
			continue;
		}

		//More than 128 bits: a prime that large cannot be a factor below 2^32.
		failed = num_pending + 1 > MAX_PENDING || mpz_probab_prime_p(c, 1) != 0 || !rho_mpz(c, d);
		if (!failed)
		{
			mpz_divexact(c, c, d);
			mpz_init_set(pending[num_pending++], d);
		}
	}

	while (num_pending > 0)
	{
		mpz_clear(pending[--num_pending]);
	}
	mpz_clears(n, d, NULL);

	if (failed)
	{
		return -1;
	}
	std::sort(factors, factors + count);
	return count;
}

//Find a nontrivial factor of an odd composite n < 2^64 with Pollard-Brent rho. Returns false if none was found.
bool split_composite64(unsigned long long n, unsigned long long* factor)
{
	return rho64(n, factor);
}

//Strong probable prime test to base 2 for odd n, 2^32 <= n < 2^64. one is 2^64 mod n, the Montgomery form of 1.
static bool sprp2_64(u64 n)
{
	u64 ninv = montgomery_inverse(n);
	u64 one = (0 - n) % n;
	u64 minus_one = n - one;

	unsigned long s;
	_BitScanForward64(&s, n - 1);
	u64 d = (n - 1) >> s;
	unsigned long top;
	_BitScanReverse64(&top, d);

	//2^d mod n, left to right. Multiplying by the base 2 is a doubling.
	u64 x = one;
	for (int i = (int)top; i >= 0; i--)
	{
		x = montgomery_multiply64(x, x, n, ninv);

		//Double x mod n if the bit is set, without branches.
		u64 doubled, reduced;
		unsigned char carry = _addcarry_u64(0, x, x, &doubled);
		unsigned char borrow = _subborrow_u64(0, doubled, n, &reduced);
		u64 keep = 0 - (u64)(borrow & (carry == 0)); //all ones if 2x < n.
		reduced ^= (reduced ^ doubled) & keep;
		u64 take = 0 - ((d >> i) & 1);
		x ^= (x ^ reduced) & take;
	}

	if (x == one || x == minus_one)
	{
		return true;
	}
	for (unsigned long i = 1; i < s; i++)
	{
		x = montgomery_multiply64(x, x, n, ninv);
		if (x == minus_one)
		{
			return true;
		}
	}
	return false;
}

//x = 2 * x mod n, for x < n.
static inline u128 double_mod(const u128& x, const u128& n)
{
	u128 r;
	unsigned char carry = _addcarry_u64(0, x.lo, x.lo, &r.lo);
	carry = _addcarry_u64(carry, x.hi, x.hi, &r.hi);
	return reduce_once(r, carry, n);
}

//Strong probable prime test to base 2 for odd n, 2^64 <= n < 2^128.
static bool sprp2_128(const u128& n)
{
	u64 ninv = montgomery_inverse(n.lo);

	//one = 2^128 mod n, the Montgomery form of 1, by doubling 2^(bits - 1) < n up to 2^128.
	unsigned long top;
	_BitScanReverse64(&top, n.hi);
	u128 one = { 0, (u64)1 << top };
	for (unsigned long i = top + 64; i < 128; i++)
	{
		one = double_mod(one, n);
	}
	u128 minus_one = subtract(n, one);

	u128 d = { n.lo - 1, n.hi }; //n is odd, so no borrow.
	unsigned long s;
	if (d.lo != 0)
	{
		_BitScanForward64(&s, d.lo);
	}
	else
	{
		_BitScanForward64(&s, d.hi);
		s += 64;
	}
	d = shift_right(d, s);

	int d_top;
	if (d.hi != 0)
	{
		_BitScanReverse64(&top, d.hi);
		d_top = (int)top + 64;
	}
	else
	{
		_BitScanReverse64(&top, d.lo);
		d_top = (int)top;
	}

	u128 x = one;
	for (int i = d_top; i >= 0; i--)
	{
		x = montgomery_square(x, n, ninv);
		u64 bit = (i >= 64) ? (d.hi >> (i - 64)) : (d.lo >> i);
		u64 take = 0 - (bit & 1);
		u128 doubled = double_mod(x, n);
		x.lo ^= (x.lo ^ doubled.lo) & take;
		x.hi ^= (x.hi ^ doubled.hi) & take;
	}

	if (equal(x, one) || equal(x, minus_one))
	{
		return true;
	}
	for (unsigned long i = 1; i < s; i++)
	{
		x = montgomery_square(x, n, ninv);
		if (equal(x, minus_one))
		{
			return true;
		}
	}
	return false;
}

bool is_probable_prime(mpz_srcptr value)
{
	int size = value->_mp_size;
	if (size > 2)
	{
		return mpz_probab_prime_p(value, 1) != 0;
	}
	if (size <= 0)
	{
		return false;
	}

	u128 n = { value->_mp_d[0], (size > 1) ? value->_mp_d[1] : 0 };
	if ((n.lo & 1) == 0)
	{
		return n.hi == 0 && n.lo == 2;
	}
	if (n.hi != 0)
	{
		return sprp2_128(n);
	}
	if (n.lo < ((u64)1 << 32))
	{
		return is_prime32(n.lo);
	}
	return sprp2_64(n.lo);
}
