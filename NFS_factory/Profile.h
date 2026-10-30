#pragma once
//Simple timing of the main parts of the program, to see where the time is spent.
//Uses the processor time stamp counter, which costs only a few cycles to read.

#include <intrin.h>

enum profile_phase {
	PROF_SETUP,           //prime product and memory setup
	PROF_READ,            //reading and parsing the relation file
	PROF_POLY_EVAL,       //polynomial values of the relations
	PROF_PRODUCT_TREE,    //product tree of the polynomial values
	PROF_ROOT,            //root of the remainder tree, (P mod root) * 2^k / root
	PROF_REMAINDER_TREE,  //scaled remainder tree down to the parents of the leaves
	PROF_LEAVES,          //remainders of the leaves, squarings and gcd
	PROF_COFACTOR,        //checking the cofactors for large primes (rho, tinyqs)
	PROF_FACTOR,          //factoring the smooth relations for output
	PROF_OTHER_SIDE,      //rebuilding the factor list of the other side (format: cado)
	PROF_WRITE,           //formatting and writing the relations
	PROF_3LP,             //the 3LP second stage, all of its phases
	PROF_NUM_PHASES
};

extern unsigned long long profile_cycles[PROF_NUM_PHASES];
extern unsigned long long profile_cycles_3lp[PROF_NUM_PHASES]; //the phases of the 3LP second stage, which are counted as one phase in profile_cycles.

inline unsigned long long profile_now()
{
	return __rdtsc();
}

//Add the time since start to a phase, and return the current time.
inline unsigned long long profile_add(profile_phase phase, unsigned long long start)
{
	unsigned long long now = __rdtsc();
	profile_cycles[phase] += now - start;
	return now;
}

void print_profile(unsigned long long total_cycles);
