#pragma once
#include "mpir.h"

//Factor a relation value known to have only prime factors below 2^32, using native 64-bit arithmetic.
//The factors are stored in increasing order, with multiplicity. Returns the number of factors, or -1 if
//the value could not be factored this way (more than 128 bits, or a factor that was not found), in which
//case the caller should use a general method.
int factor_smooth(unsigned int* factors, mpz_srcptr value, int max_factors);

//Factor a value of any size whose prime factors are all below 2^32. Up to 128 bits this is factor_smooth. Above that, or if
//factor_smooth fails, the primes below 2^16 are found by trial division, and further factors are split off with Pollard rho
//in MPIR arithmetic until the parts fit in 128 bits. The factors are stored in
//increasing order, with multiplicity. Returns the number of factors (0 for 1), or -1 if the value is 0, has a prime factor
//of 2^32 or more, has more than max_factors factors, or could not be factored.
int factor_general(unsigned int* factors, mpz_srcptr value, int max_factors);

//Find a nontrivial factor of an odd composite n < 2^64. Returns false if none was found.
bool split_composite64(unsigned long long n, unsigned long long* factor);

//Primality test for a positive number. Below 2^32 the answer is exact. Up to 128 bits it is a strong probable prime test
//to base 2 with native arithmetic, which never calls a prime composite, and calls very few composites prime. Larger numbers
//use MPIR.
bool is_probable_prime(mpz_srcptr value);
