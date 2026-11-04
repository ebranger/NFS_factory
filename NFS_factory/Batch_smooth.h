#pragma once
#include "mpirxx.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <deque>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <functional>
#include "Polynomial.h"
#include "Tree_backend.h"
#include "Profile.h"
using namespace std;

#define MAX_SAVED_PRIMES 64 //Most primes in the saved factor list of a relation from a binary relation file.
#define MAX_COFACTOR_PRIMES 96 //Most prime factors of a 3LP cofactor, which is less than 2^96.

//The relations of one batch, from reading them to writing the smooth ones: the input, the tree with their values as leaves,
//and the results of the checks.
struct Batch_slot
{
	Tree_backend* tree; //the product and remainder trees, with the relation values as leaves.
	long long *a_vals, *b_vals;
	char** list_of_factors; //holds the already known smooth side factor list, as text or (binary_lists) as the number of primes followed by the primes.
	bool binary_lists;
	bool* is_smooth; //set by the batch check for relations that are smooth.
	unsigned long long *lp1, *lp2; //the large primes of the relations
	int entries;
	long long sequence; //the number of the batch
	long long position[2]; //the read position of its last relation, as given by Set_read_position
	size_t root_bits; //of the product tree
	bool tree_factors; //the tree factored the smooth relations of the batch (Tree_backend::Leaf_factors)
	bool rebuilds_other_side; //the factor list of the other side is rebuilt for some of the relations (see Batch_smooth)
	mpz_class temp_number;
	vector<unsigned char> leaf_result; //cofactor_result of each leaf
	vector<string> relation_texts;
	vector<unsigned char> relation_result;
	vector<int> candidate_index; //3LP candidates of the batch, sent on after the batch is checked.
	vector<mpz_class> candidate_cofactor;
};

class Batch_smooth
{
public:
	//other_side_poly: the polynomial of the side that is not batch checked. The factor list of that side is rebuilt from its norm, with
	//all primes and their multiplicities, for every relation if rebuild_other_side (as CADO-NFS requires), and otherwise for the
	//relations that were added without one (Add_number_binary with no primes, from a file of only a and b); the others are printed as
	//they were read, if they divide the norm. The primes below print_limit are left out. On the GPU, the norms are factored with the
	//relation values.
	//backend: the implementation of the tree computations.
	//three_lp_stage: if true, this is the 3LP second stage. Its inputs are the cofactors given with Add_candidate, which are checked against all primes below 2^lbp (fblim and mfb are not used).
	//pipeline: the number of batches in progress at once (see Batch_slot), each with its own tree, or 0 for the default of the
	//backend (default_pipeline). 1 checks each batch when it is full, in the calling thread. With more, the batches are checked in a pipeline of threads while the next ones are read: the
	//polynomial values and the trees (by up to four threads, as the trees use the GPU), and then the cofactors and the output (one
	//thread, in the order of the batches, so that the output does not depend on the pipeline). The second stage always uses 1.
	Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side, int print_limit, Polynomial* other_side_poly, bool rebuild_other_side, int max_tree_levels, tree_backend_kind backend, bool three_lp_stage = false, int pipeline = 1);
	~Batch_smooth(void);
	void Finish(); //Check the remaining relations, and any 3LP candidates in the second stage. Call once, after the last relation.
	bool Add_number(long long a_val, long long b_val, const char* relation_factor_list);
	bool Add_number_binary(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes);
	bool Add_candidate(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes, mpz_srcptr cofactor);
	long long Get_num_relations_found(); //including the 3LP second stage
	long long Get_num_3lp_relations_found(); //by the 3LP second stage
	long long Get_num_relations_checked(); //relations of the batches checked so far (not counting the second stage)
	int Get_batch_size() { return num_relations_per_batch; }
	int Get_tree_levels() { return tree_size; }
	int Get_pipeline() { return (int)slots.size(); }

	//The read position of the next relations added (two numbers for the caller, such as file and record), which is given back
	//by the batch done callback.
	void Set_read_position(long long first, long long second);
	//Called after each batch, in their order, when its relations are written to the output file (in the thread that writes them),
	//with the read position of its last relation, and the size of the output file then (-1 if it could not be written).
	void Set_batch_done(function<void(const long long position[2], long long output_size)> callback);
	//The time the thread that adds the relations spent in Add_* and Finish waiting for the checks, or doing them (pipeline 1).
	unsigned long long Get_reader_blocked_cycles();

	//3LP candidates, cofactors from 2 lbp + 1 to mfb bits, are passed to a second stage that this object then owns (NULL: not used),
	//and/or appended to a candidate file. Returns false if the file could not be opened.
	bool Set_three_lp_handling(Batch_smooth* second_stage, string candidate_filename);

	int side;

private:
	mpz_class* prime_product;
	void Get_saved_list(Batch_slot& s, int index, string& list);
	bool Saved_list_missing(Batch_slot& s, int index);
	void Get_saved_primes(Batch_slot& s, int index, vector<unsigned long long>& primes);

	void Setup_prime_product(long long smooth_bound, bool use_gpu);
	Batch_slot* New_slot(int tree_levels, int max_bits_in_input, tree_backend_kind backend);
	void Delete_slot(Batch_slot* s);

	//The steps of a batch: the trees, then the checks of the cofactors and the output.
	void Submit();
	void Run_trees(Batch_slot& s);
	void Run_checks(Batch_slot& s);
	void Tree_worker();
	void Check_worker();
	void Output_solutions(Batch_slot& s);
	enum relation_status { RELATION_WRITTEN, RELATION_FACTORING_FAILED, RELATION_OTHER_SIDE_FAILED };
	relation_status Format_relation(Batch_slot& s, int index, string& relation_text);
	bool Complete_other_side(Batch_slot& s, int index, string& list);
	bool Saved_list_divides(Batch_slot& s, int index);

	string save_filename;
	FILE* output; //open from the first batch on; flushed after each batch, so that its size is that of the relations written
	int tree_size;
	int max_entries_before_check;
	atomic<long long> num_relations_found;
	atomic<long long> num_relations_checked;
	long long num_factoring_failures;
	long long num_other_side_failures;
	int num_relations_per_batch;
	int min_saved_prime_limit;

	Polynomial* poly;
	Polynomial* other_poly; //for rebuilding the other side's factor list
	bool rebuild_all; //for every relation, not only those without a saved list
	bool other_side_on_tree; //the tree factors the other side (Tree_backend::Set_other_side)

	bool using_snfs_deg4_binomial;
	bool binomial_fits; //the coefficients of the binomial fit in 63 bits, for Binomial_value_native.
	long long binomial_c4, binomial_c0;
	bool Binomial_value_native(mpz_ptr result, long long a, long long b);
	long long max_batch_prime;

	int one_lp_bit_limit;
	long long two_lp_min_limit;
	int max_cofactor_bit_size;
	bool do_cofactorization;
	//The cofactor checks of the leaves: Classify_cofactor is thread safe and does all but tinyqs, which is not (it keeps state
	//between calls), so the cofactors that need it are split afterwards in leaf order, as before.
	enum cofactor_result { COFACTOR_NOT_SMOOTH, COFACTOR_SMOOTH, COFACTOR_CANDIDATE, COFACTOR_TINYQS };
	cofactor_result Classify_cofactor(Batch_slot& s, mpz_class& number, unsigned int index);
	cofactor_result Check_two_primes(Batch_slot& s, unsigned long long p, unsigned long long q, unsigned int index);
	bool Split_with_tinyqs(Batch_slot& s, mpz_class& number, unsigned int index);
	mpz_class *factor1, *factor2;

	//3LP handling.
	bool three_lp_stage; //this object is the second stage.
	int three_lp_max_bits; //largest 3LP candidate, mfb if it is more than 2 lbp, otherwise 0 (no 3LP).
	mpz_class three_lp_min; //smallest 3LP candidate, fblim^3.
	Batch_smooth* second_stage;
	FILE* candidate_file;
	long long num_candidates;
	long long num_candidates_skipped; //could not be saved (a, b or saved primes too large for the candidate format).
	void Send_candidates(Batch_slot& s);

	//The pipeline. A slot is free, being filled (filling), waiting for or in the tree stage (tree_queue), or waiting for or in
	//the checks, which take the batches in sequence.
	profile_counter* profile; //profile_cycles, or profile_cycles_3lp for the second stage
	vector<Batch_slot*> slots;
	Batch_slot* filling;
	long long read_position[2];
	function<void(const long long position[2], long long output_size)> batch_done;
	long long output_size;
	bool pipelined;
	mutex lock;
	condition_variable changed;
	deque<Batch_slot*> free_slots;
	deque<Batch_slot*> tree_queue;
	vector<Batch_slot*> trees_done; //waiting for the checks
	long long submitted; //batches submitted
	long long checked; //batches checked and written
	bool stopping;
	vector<thread> workers;
	unsigned long long reader_blocked;
};
