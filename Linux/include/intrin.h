#pragma once
//Linux stand-in for MSVC's <intrin.h>: the intrinsics NFS_factory uses, for GCC and Clang on x86-64.

#include <x86intrin.h>

//64 x 64 -> 128 bit product: returns the low half, and the high half in *high. uint64_t is unsigned long on Linux, so both
//pointer types are accepted.
static inline unsigned long long _umul128(unsigned long long a, unsigned long long b, unsigned long long* high)
{
	unsigned __int128 p = (unsigned __int128)a * b;
	*high = (unsigned long long)(p >> 64);
	return (unsigned long long)p;
}

static inline unsigned long long _umul128(unsigned long long a, unsigned long long b, unsigned long* high)
{
	unsigned __int128 p = (unsigned __int128)a * b;
	*high = (unsigned long)(p >> 64);
	return (unsigned long long)p;
}

//GCC's _addcarry_u64 takes an unsigned long long pointer; this one takes uint64_t's unsigned long.
static inline unsigned char _addcarry_u64(unsigned char carry, unsigned long long a, unsigned long long b, unsigned long* sum)
{
	unsigned long long s;
	carry = _addcarry_u64(carry, a, b, &s);
	*sum = (unsigned long)s;
	return carry;
}

//The index of the lowest (highest) set bit of mask; return 0 if mask is zero, and leave *index alone.
static inline unsigned char _BitScanForward64(unsigned long* index, unsigned long long mask)
{
	if (mask == 0)
	{
		return 0;
	}
	*index = (unsigned long)__builtin_ctzll(mask);
	return 1;
}

static inline unsigned char _BitScanReverse64(unsigned long* index, unsigned long long mask)
{
	if (mask == 0)
	{
		return 0;
	}
	*index = 63 - (unsigned long)__builtin_clzll(mask);
	return 1;
}

//128 / 64 -> 64 bit division of high:low by divisor, which must be more than high: returns the quotient, and the remainder in
//*remainder.
static inline unsigned long long _udiv128(unsigned long long high, unsigned long long low, unsigned long long divisor, unsigned long long* remainder)
{
	unsigned long long quotient, rest;
	__asm__("divq %4" : "=a"(quotient), "=d"(rest) : "a"(low), "d"(high), "rm"(divisor));
	*remainder = rest;
	return quotient;
}

static inline unsigned long long _udiv128(unsigned long long high, unsigned long long low, unsigned long long divisor, unsigned long* remainder)
{
	unsigned long long rest;
	unsigned long long quotient = _udiv128(high, low, divisor, &rest);
	*remainder = (unsigned long)rest;
	return quotient;
}
