#pragma once
#include "mpirxx.h"
#include "Polynomial.h"
using namespace std;

void setup_prime_product(Polynomial* poly, uint32 max_prime, mpz_class* prime_product);
