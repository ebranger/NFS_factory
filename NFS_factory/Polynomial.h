#pragma once

#include "mpirxx.h"
#include <stdint.h>

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef uint64_t uint64;
	
typedef char int8;
typedef short int16;
typedef int int32;
typedef int64_t int64;

class Polynomial
{
public:
	Polynomial();
	Polynomial(Polynomial* input);
	~Polynomial(void);

	void copy(Polynomial* input);

	uint8 get_degree();

	void set_coeff(uint8 coeff, mpz_class* value);
	void set_coeff(uint8 coeff, int32 value);

	void clear();

	void evaluate_polynomial(mpz_class* result, mpz_class* value);
	void evaluate_polynomial(mpz_class* result, uint32 value);

	void evaluate_polynomial_modulo(mpz_class* result, mpz_class* value, mpz_class* module);
	void evaluate_polynomial_modulo(mpz_class* result, uint32 value, uint32 module);

	void evaluate_polynomial_homogenous(mpz_class* result, int32 a, int32 b);
	void evaluate_polynomial_homogenous(mpz_class* result, int64 a, int64 b);

	bool poly_is_binomial_SNFS_deg4();
	void evaluate_binomial_polynomial_homogeneous_SNFS_deg4(mpz_class* result, int64 a, int64 b, mpz_class* scratchpad_mpz);

	void debug_print_poly();
    void update_degree();
	void mult_modpoly(Polynomial* a, Polynomial* modpoly, uint32 p);
	void reduce_modp(uint32 p);
	void reduce_modp(mpz_class* p);
	void make_monic(uint32 p);
	void make_monic(mpz_class* prime);
	void reduce_modpoly(Polynomial* modulus, uint32 p);
	void gcd(Polynomial* a, Polynomial* b, uint32 p);
	void square_poly_modpoly(Polynomial* modpoly, uint32 p);
	void poly_exponentiate(uint32 power, Polynomial* modpoly, uint32 p);
	void divide_out_root(uint32 root, uint32 p);
	void get_polyroots_modp(uint32* numroots, uint32* rootlist, uint32 p);
private:
	uint8 degree;
	mpz_class** coefficients;
};

