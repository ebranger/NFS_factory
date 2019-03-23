#pragma once

#define s32 int
#define u32 unsigned int
#define s64 long long int //modified to be long long, i.e.64 bit.
#define u64 unsigned long long int
#define s16 short int
#define u16 unsigned short int
#define u8 unsigned char

#include "mpir.h"

u32 squfof(mpz_t n);
