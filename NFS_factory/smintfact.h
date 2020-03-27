#pragma once

#define s32 int
#define u32 unsigned int
extern "C" {
	int prho(mpz_t factor1, mpz_t factor2, mpz_t n, s32 c, s32 maxIts);
	int factor(u32 *factors, mpz_t n, int useTrialDivision);
	u32 squfof(mpz_t n);
}
