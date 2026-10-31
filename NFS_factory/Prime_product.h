#pragma once
#include "mpirxx.h"
#include "Polynomial.h"
using namespace std;

//Set up the product of all primes up to max_prime that can divide the homogeneous norm F(a,b) of the polynomial for coprime a, b.
//Small primes are included to higher powers, to find relations where they occur to a higher power.
//For a polynomial of degree 1 every prime can divide the norm, so all primes are used. The primes are checked in parallel (with
//the threads option). With use_gpu (in the GPU build), the product is computed on the GPU; it is the same number.
void setup_prime_product(Polynomial* poly, uint32 max_prime, mpz_class* prime_product, bool use_gpu);
