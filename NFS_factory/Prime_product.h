#pragma once
#include "mpirxx.h"
#include "Polynomial.h"
using namespace std;

//Set up the product of all primes up to max_prime that can divide the homogeneous norm F(a,b) of the polynomial for coprime a, b.
//Small primes are included to higher powers, to find relations where they occur to a higher power.
//For a polynomial of degree 1 every prime can divide the norm, so all primes are used.
void setup_prime_product(Polynomial* poly, uint32 max_prime, mpz_class* prime_product);
