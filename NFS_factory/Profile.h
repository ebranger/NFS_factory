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

//The time that the batches of a pipeline spend in each of its stages, added up over the batches, for the batches of the relations
//(index 0) and those of the batch check of the other side's norms (index 1), with the number of batches: what a batch waits for.
enum profile_stage {
	STAGE_WAIT_TREES,  //read, waiting for a thread of the trees
	STAGE_TREES,       //in the trees
	STAGE_WAIT_CHECKS, //waiting for a thread of the checks
	STAGE_COFACTORS,   //waiting for the GPU to check the cofactors (rho), and their checks on the CPU
	STAGE_FACTORS,     //waiting for the GPU to factor the relations found
	STAGE_RELATIONS,   //the relations found, on the CPU
	STAGE_WAIT_ORDER,  //waiting for the batches before it
	STAGE_ORDER,       //the steps in the order of the batches, and handing the relations to the writer
	STAGE_NUM
};
extern profile_counter profile_stages[2][STAGE_NUM];
extern std::atomic<long long> profile_stage_batches[2];

//seconds: the time of the run, for the times in ms (0: not printed).
void print_profile(unsigned long long total_cycles, double seconds = 0);

//The output during a run: a status line, which on a console is rewritten in place (at most twice a second), and otherwise written
//as a line of its own at most once a minute, and messages, which are printed above the status line. Thread safe.
extern bool verbose_output; //also print a line for each batch (-verbose)
bool console_is_terminal();
void console_status(const std::string& text);
void console_message(const std::string& text);
//The last status line, always written, after which messages are printed as usual.
void console_end_status(const std::string& text);
