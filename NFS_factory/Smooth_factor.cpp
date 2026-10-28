//Factoring of smooth relation values with native 64-bit arithmetic.
//After the large primes are divided out, the value of a relation from the batch check only has prime factors below
//the factor base limit, which is less than 2^32. The value then fits in 128 bits, and any remaining cofactor of at
//least 2^32 must be composite. Small primes are found by trial division, the rest by Pollard-Brent rho with Montgomery
//multiplication modulo the (up to 128-bit) cofactor.

#include "Smooth_factor.h"
#include <intrin.h>
#include <algorithm>

typedef unsigned long long u64;
typedef unsigned int u32;

#define RHO_MAX_ITERATIONS (1 << 20) //sqrt(2^32) = 65536 iterations are expected for the largest primes, so this is plenty.
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
	if (t2 != 0 || !less(r, n))
	{
		r = subtract(r, n);
	}
	return r;
}

//The rho map y -> y^2 + c mod n. Working on Montgomery representations is fine, since the result mod any prime p
//dividing n only depends on y mod p.
static inline u128 rho_step(const u128& y, const u128& c, const u128& n, u64 ninv)
{
	u128 r = montgomery_multiply(y, y, n, ninv);
	u128 sum;
	unsigned char carry = _addcarry_u64(0, r.lo, c.lo, &sum.lo);
	carry = _addcarry_u64(carry, r.hi, c.hi, &sum.hi);
	if (carry || !less(sum, n))
	{
		sum = subtract(sum, n);
	}
	return sum;
}

static inline u128 difference(const u128& a, const u128& b)
{
	return less(a, b) ? subtract(b, a) : subtract(a, b);
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
		else if (!rho(c, &d) || d.hi != 0)
		{
			return -1;
		}

		u64 remainder;
		pending[num_pending++] = d;
		pending[num_pending++] = divide(c, d.lo, &remainder);
	}

	std::sort(factors, factors + count);
	return count;
}

//Find a nontrivial factor of an odd composite n < 2^64 with Pollard-Brent rho. Returns false if none was found.
bool split_composite64(unsigned long long n, unsigned long long* factor)
{
	return rho64(n, factor);
}
