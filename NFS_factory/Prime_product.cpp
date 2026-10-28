#include "Prime_product.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
using namespace std;

//Only primes that can divide the norm F(a,b) of a relation are needed in the prime product. For a polynomial of degree d,
//a prime p can divide F(a,b) with a, b coprime only if f has a root mod p, or p divides the leading coefficient.
//For a typical degree 4 polynomial that is less than half of the primes, so the prime product gets much smaller.
//All arithmetic mod p is done with 64-bit integers, since p < 2^32.

#define SIEVE_SEGMENT 32768 //Odd numbers per segment of the sieve of Eratosthenes, small enough to fit in the cache.
#define MAX_POLY_DEGREE 10
#define BRUTE_FORCE_LIMIT 100 //For primes below this, roots are found by trying all values.
#define SMALL_PRIME_POWERS "10821610800" //The smallest primes are included to a higher power.

//Find all primes up to max_prime, using a segmented sieve of Eratosthenes over the odd numbers.
static void find_primes(uint32 max_prime, vector<uint32>& primes)
{
	primes.clear();
	if (max_prime < 2)
	{
		return;
	}
	primes.push_back(2);

	//The primes up to sqrt(max_prime) are used for sieving.
	uint32 root = (uint32)sqrt((double)max_prime);
	while ((uint64)root * root > max_prime)
	{
		root--;
	}
	while ((uint64)(root + 1) * (root + 1) <= max_prime)
	{
		root++;
	}

	vector<uint8> is_prime(root + 1, 1);
	vector<uint32> sieve_primes;
	for (uint32 i = 3; i <= root; i += 2)
	{
		if (is_prime[i])
		{
			sieve_primes.push_back(i);
			for (uint64 j = (uint64)i * i; j <= root; j += 2 * i)
			{
				is_prime[(size_t)j] = 0;
			}
		}
	}

	//Each segment covers the odd numbers low, low + 2, ..., low + 2 * (SIEVE_SEGMENT - 1).
	vector<uint8> segment(SIEVE_SEGMENT);
	vector<uint64> next_multiple(sieve_primes.size()); //Next odd multiple of each sieving prime to cross out.
	for (size_t i = 0; i < sieve_primes.size(); i++)
	{
		next_multiple[i] = (uint64)sieve_primes[i] * sieve_primes[i];
	}

	for (uint64 low = 3; low <= max_prime; low += 2 * (uint64)SIEVE_SEGMENT)
	{
		uint64 high = low + 2 * (uint64)SIEVE_SEGMENT;
		fill(segment.begin(), segment.end(), 1);

		for (size_t i = 0; i < sieve_primes.size(); i++)
		{
			uint64 step = 2 * (uint64)sieve_primes[i];
			uint64 j = next_multiple[i];
			for (; j < high; j += step)
			{
				segment[(size_t)((j - low) / 2)] = 0;
			}
			next_multiple[i] = j;
		}

		for (uint32 k = 0; k < SIEVE_SEGMENT; k++)
		{
			uint64 n = low + 2 * (uint64)k;
			if (n > max_prime)
			{
				break;
			}
			if (segment[k])
			{
				primes.push_back((uint32)n);
			}
		}
	}
}

static inline uint32 mulmod(uint32 a, uint32 b, uint32 p)
{
	return (uint32)((uint64)a * b % p);
}

static uint32 powmod(uint32 base, uint32 exponent, uint32 p)
{
	uint64 result = 1;
	uint64 b = base % p;

	while (exponent > 0)
	{
		if (exponent & 1)
		{
			result = result * b % p;
		}
		b = b * b % p;
		exponent >>= 1;
	}
	return (uint32)result;
}

static uint32 gcd(uint32 a, uint32 b)
{
	while (b != 0)
	{
		uint32 t = a % b;
		a = b;
		b = t;
	}
	return a;
}

//Degree of a polynomial, or -1 for the zero polynomial.
static int poly_degree(const uint32* a, int max_degree)
{
	for (int i = max_degree; i >= 0; i--)
	{
		if (a[i] != 0)
		{
			return i;
		}
	}
	return -1;
}

//result = a * b mod (g, p), where g is monic of degree d and a, b have degree less than d.
static void poly_mulmod(const uint32* a, const uint32* b, const uint32* g, int d, uint32 p, uint32* result)
{
	uint64 product[2 * MAX_POLY_DEGREE];

	for (int i = 0; i < 2 * d - 1; i++)
	{
		product[i] = 0;
	}

	for (int i = 0; i < d; i++)
	{
		for (int j = 0; j < d; j++)
		{
			product[i + j] = (product[i + j] + (uint64)a[i] * b[j]) % p;
		}
	}

	//Reduce using x^d = -(g[d-1] x^(d-1) + ... + g[0])
	for (int k = 2 * d - 2; k >= d; k--)
	{
		uint64 c = product[k];
		if (c != 0)
		{
			for (int j = 0; j < d; j++)
			{
				product[k - d + j] = (product[k - d + j] + (p - c) * (uint64)g[j]) % p;
			}
		}
	}

	for (int i = 0; i < d; i++)
	{
		result[i] = (uint32)product[i];
	}
}

//Degree of gcd(a, b) mod p. a and b are overwritten.
static int poly_gcd_degree(uint32* a, uint32* b, int max_degree, uint32 p)
{
	int da = poly_degree(a, max_degree);
	int db = poly_degree(b, max_degree);

	while (db >= 0)
	{
		//a = a mod b
		uint32 inverse = powmod(b[db], p - 2, p);
		while (da >= db)
		{
			uint32 c = mulmod(a[da], inverse, p);
			for (int j = 0; j <= db; j++)
			{
				a[da - db + j] = (uint32)((a[da - db + j] + (uint64)(p - c) * b[j]) % p);
			}
			da = poly_degree(a, da - 1);
		}

		swap(a, b);
		swap(da, db);
	}
	return da;
}

//Check if the polynomial with coefficients f[0..d] (mod p, f[d] nonzero) has a root mod p.
static bool has_root(const uint32* f, int d, uint32 p)
{
	if (f[0] == 0)
	{
		//x = 0 is a root.
		return true;
	}

	if (p < BRUTE_FORCE_LIMIT)
	{
		for (uint32 x = 1; x < p; x++)
		{
			uint64 value = 0;
			for (int i = d; i >= 0; i--)
			{
				value = (value * x + f[i]) % p;
			}
			if (value == 0)
			{
				return true;
			}
		}
		return false;
	}

	//Make the polynomial monic.
	uint32 g[MAX_POLY_DEGREE + 1];
	uint32 inverse = powmod(f[d], p - 2, p);
	for (int i = 0; i <= d; i++)
	{
		g[i] = mulmod(f[i], inverse, p);
	}

	//Compute x^p mod g by squaring and multiplying, starting from x.
	uint32 r[MAX_POLY_DEGREE + 1] = { 0 };
	uint32 x[MAX_POLY_DEGREE + 1] = { 0 };
	uint32 temp[MAX_POLY_DEGREE + 1];
	r[1] = 1;
	x[1] = 1;

	int top_bit = 31;
	while (((p >> top_bit) & 1) == 0)
	{
		top_bit--;
	}

	for (int bit = top_bit - 1; bit >= 0; bit--)
	{
		poly_mulmod(r, r, g, d, p, temp);
		if ((p >> bit) & 1)
		{
			poly_mulmod(temp, x, g, d, p, r);
		}
		else
		{
			copy(temp, temp + d, r);
		}
	}

	//f has a root mod p if gcd(g, x^p - x) is not constant.
	r[1] = (r[1] + p - 1) % p;
	for (int i = d; i <= MAX_POLY_DEGREE; i++)
	{
		r[i] = 0;
	}

	if (poly_degree(r, d - 1) < 0)
	{
		//x^p - x is divisible by g, so every root of g is in the field.
		return true;
	}

	return poly_gcd_degree(g, r, d, p) >= 1;
}

//Check if a binomial f = c_d x^d + c_0 has a root mod p, i.e. if -c_0 / c_d is a d-th power mod p.
static bool binomial_has_root(uint32 c_d, uint32 c_0, int d, uint32 p)
{
	if (c_0 == 0)
	{
		return true;
	}

	uint32 k = mulmod(p - c_0, powmod(c_d, p - 2, p), p);
	uint32 g = gcd((uint32)d, p - 1);
	return powmod(k, (p - 1) / g, p) == 1;
}

//Multiply the factors into result, using a balanced product tree.
static void multiply_list(const uint32* factors, size_t n, mpz_class& result)
{
	if (n <= 64)
	{
		result = 1;
		uint64 accumulated = 1;

		for (size_t i = 0; i < n; i++)
		{
			if (accumulated > 0xFFFFFFFFULL)
			{
				mpz_mul_ui(result.get_mpz_t(), result.get_mpz_t(), accumulated);
				accumulated = 1;
			}
			accumulated *= factors[i];
		}
		mpz_mul_ui(result.get_mpz_t(), result.get_mpz_t(), accumulated);
		return;
	}

	mpz_class left, right;
	multiply_list(factors, n / 2, left);
	multiply_list(factors + n / 2, n - n / 2, right);
	result = left * right;
}

void setup_prime_product(Polynomial* poly, uint32 maxprime, mpz_class* prime_product)
{
	uint32 square_limit = uint32(sqrt((double)maxprime));
	uint32 fourth_power_limit = uint32(sqrt(sqrt((double)maxprime)));

	mpz_class small_primes;
	small_primes = SMALL_PRIME_POWERS;

	if (poly->get_degree() <= 1)
	{
		//Degree 1 polynomial, all primes can divide the norm. Just do a primorial.
		mpz_class temp_number;
		mpz_primorial_ui(prime_product->get_mpz_t(), maxprime);
		mpz_primorial_ui(temp_number.get_mpz_t(), square_limit); //Extra small primes. this saves having to do some squarings in the batch tree bottom node, and adds very few bits to the product.
		*prime_product = *prime_product * temp_number;
		mpz_primorial_ui(temp_number.get_mpz_t(), fourth_power_limit); // Extra even smaller primes.
		*prime_product = *prime_product * (temp_number * temp_number);
		*prime_product = *prime_product * (small_primes * small_primes); //and the smallest primes to an even higher level!
		return;
	}

	int d = poly->get_degree();

	//Check if the coefficients fit in 64 bits, then they can be reduced mod p faster.
	bool small_coefficients = true;
	long long coefficients[MAX_POLY_DEGREE + 1];
	for (int i = 0; i <= d; i++)
	{
		if (mpz_fits_si_p(poly->get_coefficient(i).get_mpz_t()))
		{
			coefficients[i] = mpz_get_si(poly->get_coefficient(i).get_mpz_t());
		}
		else
		{
			small_coefficients = false;
		}
	}

	bool binomial = true;
	for (int i = 1; i < d; i++)
	{
		if (poly->get_coefficient(i) != 0)
		{
			binomial = false;
		}
	}

	vector<uint32> primes;
	find_primes(maxprime, primes);

	//The primes that can divide the norm, with small primes repeated to include them to higher powers.
	vector<uint32> factors;
	size_t num_dividing_primes = 0;
	uint32 f[MAX_POLY_DEGREE + 1];

	for (size_t n = 0; n < primes.size(); n++)
	{
		uint32 p = primes[n];

		for (int i = 0; i <= d; i++)
		{
			if (small_coefficients)
			{
				long long c = coefficients[i] % (long long)p;
				f[i] = (uint32)(c < 0 ? c + p : c);
			}
			else
			{
				f[i] = (uint32)mpz_fdiv_ui(poly->get_coefficient(i).get_mpz_t(), p);
			}
		}

		bool can_divide;
		if (f[d] == 0)
		{
			//p divides the leading coefficient, so it divides F(a,b) when it divides b.
			can_divide = true;
		}
		else if (binomial && p >= BRUTE_FORCE_LIMIT)
		{
			can_divide = binomial_has_root(f[d], f[0], d, p);
		}
		else
		{
			can_divide = has_root(f, d, p);
		}

		if (can_divide)
		{
			num_dividing_primes++;
			factors.push_back(p);
			if (p <= square_limit)
			{
				factors.push_back(p);
			}
			if (p <= fourth_power_limit)
			{
				factors.push_back(p);
				factors.push_back(p);
			}
		}
	}

	multiply_list(factors.data(), factors.size(), *prime_product);
	*prime_product = *prime_product * (small_primes * small_primes);

	cout << num_dividing_primes << " of the " << primes.size() << " primes up to " << maxprime << " can divide the norm." << endl;
}
