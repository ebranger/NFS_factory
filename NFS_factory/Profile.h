#pragma once
//Simple timing of the main parts of the program, to see where the time is spent.
//Uses the processor time stamp counter, which costs only a few cycles to read. The counters can be added to from several threads;
//when the batches are checked in a pipeline (see Batch_smooth), the phases overlap, and each counts the time spent in it.

#include <intrin.h>
#include <atomic>
#include <string>

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
	PROF_OTHER_SIDE,      //the writer: the batch check of the other side's norms of the relations found, and finishing their lines
	PROF_WRITE,           //formatting and writing the relations
	PROF_3LP,             //the 3LP second stage, all of its phases
	//The waits of the pipeline, which are not part of the total:
	PROF_WAIT_READER,     //the reader waiting for a free batch slot
	PROF_WAIT_TREES,      //the tree stage waiting for a batch to be read
	PROF_WAIT_CHECKS,     //the cofactor and output stage waiting for the trees of the next batch
	PROF_NUM_PHASES
};
#define PROF_FIRST_WAIT PROF_WAIT_READER

typedef std::atomic<unsigned long long> profile_counter;
extern profile_counter profile_cycles[PROF_NUM_PHASES];
extern profile_counter profile_cycles_3lp[PROF_NUM_PHASES]; //the phases of the 3LP second stage, which are counted as one phase in profile_cycles.
extern profile_counter profile_cycles_other[PROF_NUM_PHASES]; //the phases of the batch check of the other side, part of PROF_OTHER_SIDE.

inline unsigned long long profile_now()
{
	return __rdtsc();
}

//Add the time since start to a phase of a set of counters, and return the current time.
inline unsigned long long profile_add(profile_counter* counters, profile_phase phase, unsigned long long start)
{
	unsigned long long now = __rdtsc();
	counters[phase] += now - start;
	return now;
}

inline unsigned long long profile_add(profile_phase phase, unsigned long long start)
{
	return profile_add(profile_cycles, phase, start);
}

void print_profile(unsigned long long total_cycles);

//The output during a run: a status line, which on a console is rewritten in place (at most twice a second), and otherwise written
//as a line of its own at most once a minute, and messages, which are printed above the status line. Thread safe.
extern bool verbose_output; //also print a line for each batch (-verbose)
bool console_is_terminal();
void console_status(const std::string& text);
void console_message(const std::string& text);
//The last status line, always written, after which messages are printed as usual.
void console_end_status(const std::string& text);
