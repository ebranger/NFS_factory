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
#include <unordered_map>
#include "Polynomial.h"
#include "Tree_backend.h"
#include "Profile.h"
using namespace std;

#define MAX_SAVED_PRIMES 64 //Most primes in the saved factor list of a relation from a binary relation file.
#define MAX_COFACTOR_PRIMES 96 //Most prime factors of a 3LP cofactor, which is less than 2^96.
#define OTHER_STAGE_BOUND (1 << 28) //The norms of the other side of the relations found are batch checked against the primes up to this
#define OTHER_STAGE_LARGE_PRIME_BITS 32 //and their cofactors split into two large primes of up to this many bits (see Batch_smooth);
#define OTHER_STAGE_MAX_COFACTOR_BITS 128 //larger cofactors of up to this many bits are factored apart from the rest of the norm.
#define OTHER_STAGE_LEVELS 18 //Levels of the tree of those norms: 2^17 of them at a time.

class Batch_smooth;

//The relations found in a batch (or a batch of the 3LP second stage), on their way to the output file: see Batch_smooth::Queue_output.
struct Output_batch
{
	Batch_smooth* source; //which counts the relations written
	vector<string> lines; //the relations, in the order of the batch; a line waiting for the factors of its other side is empty
	struct Waiting
	{
		size_t line;
		long long a, b;
		string checked_side; //the factor list of the batch checked side
		vector<unsigned long long> saved; //the saved factor list of the other side, if any, which must divide the norm
		long long ticket; //of the norm in the writer's other_stage, -1 if none
	};
	vector<Waiting> waiting;
	bool batch_end; //a batch of the object that writes: after it is written, its position goes to the batch done callback
	long long position[2];
};

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
	mpz_class temp_number;
	vector<unsigned char> leaf_result; //cofactor_result of each leaf
	vector<string> relation_texts;
	vector<unsigned char> relation_result;
	vector<unsigned char> relation_waiting; //for the factors of its other side
	vector<int> candidate_index; //3LP candidates of the batch, sent on after the batch is checked.
	vector<mpz_class> candidate_cofactor;
};

class Batch_smooth
{
public:
	//other_side_poly: the polynomial of the side that is not batch checked. The factor list of that side is rebuilt from its norm, with
	//all primes and their multiplicities, for every relation if rebuild_other_side (as CADO-NFS requires), and otherwise for the
	//relations that were added without one (Add_number_binary with no primes, from a file of only a and b); the others are printed as
	//they were read, if they divide the norm. The primes below print_limit are left out. On the GPU, the norms of the relations found are
	//batch checked in a stage of their own (see other_stage), and their lines are finished when they are written.
	//backend: the implementation of the tree computations.
	//three_lp_stage: if true, this is the 3LP second stage. Its inputs are the cofactors given with Add_candidate, which are checked against all primes below 2^lbp (fblim and mfb are not used).
	//pipeline: the number of batches in progress at once (see Batch_slot), each with its own tree, or 0 for the default of the
	//backend (default_pipeline). 1 checks each batch when it is full, in the calling thread. With more, the batches are checked in a pipeline of threads while the next ones are read: the
	//polynomial values and the trees (by up to four threads, as the trees use the GPU), then the cofactors and the factor lists of the
	//relations found (one thread), and the output (one thread), both in the order of the batches, so that the output does not depend
	//on the pipeline. The second stage always uses 1, and its relations are written by the writer of the first stage.
	Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side, int print_limit, Polynomial* other_side_poly, bool rebuild_other_side, int max_tree_levels, tree_backend_kind backend, bool three_lp_stage = false, int pipeline = 1);
	~Batch_smooth(void);
	void Finish(); //Check the remaining relations, and any 3LP candidates in the second stage. Call once, after the last relation.
	bool Add_number(long long a_val, long long b_val, const char* relation_factor_list);
	bool Add_number_binary(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes);
	bool Add_candidate(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes, mpz_srcptr cofactor);
	long long Get_num_relations_found(); //including the 3LP second stage, and those not yet written (see Set_batch_done)
	long long Get_num_3lp_relations_found(); //by the 3LP second stage
	long long Get_num_relations_checked(); //relations of the batches checked so far (not counting the second stage)
	int Get_batch_size() { return num_relations_per_batch; }
	int Get_tree_levels() { return tree_size; }
	int Get_pipeline() { return (int)slots.size(); }

	//The read position of the next relations added (two numbers for the caller, such as file and record), which is given back
	//by the batch done callback.
	void Set_read_position(long long first, long long second);
	//Called after each batch, in their order, when its relations are written to the output file (in the thread that writes them),
	//with the read position of its last relation, the size of the output file then (-1 if it could not be written), and the number
	//of relations written so far. The relations are found before they are written, by up to a tree of the batch check of the other
	//side when that is used (2^17 relations), so this can be less than Get_num_relations_found.
	void Set_batch_done(function<void(const long long position[2], long long output_size, long long written)> callback);
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
	relation_status Format_relation(Batch_slot& s, int index, string& relation_text, bool& waiting);
	relation_status Factor_checked_side(Batch_slot& s, int index, unsigned int* large, int& num_large, unsigned int* factors, int& num_factors);
	string Relation_line(long long a, long long b, const string& checked_side, const string& other_side);
	bool Complete_other_side(long long a, long long b, vector<unsigned long long>& primes, const unsigned int* factors, int num_factors, string& list);
	bool Saved_list_divides(Batch_slot& s, int index);

	//The output. The relations found in each batch go, in order, to the writer: with a pipeline its own thread, otherwise the thread
	//that checks the batch. The second stage gives its relations to the writer of the first stage (writer), so that they are written
	//between the right batches, as when the checks wrote them. The writer gives the norms of the other side that are to be factored
	//to other_stage, keeps the batches in held until their factors are there, finishes their lines, writes them, and after a batch
	//calls the batch done callback.
	Batch_smooth* writer; //this, or the first stage for the second stage
	Other_side_factoring* other_factoring; //NULL if the other side's norms are factored on the CPU; on the GPU, for those other_stage leaves
	tree_backend_kind backend_kind;
	void Queue_output(Output_batch* o);
	void Writer_take(Output_batch* o);
	void Write_held(bool flush);
	bool Batch_ready(const Output_batch& o);
	void Write_batch(Output_batch& o);
	void Writer_worker();
	void Flush_writer();
	deque<Output_batch*> output_queue;
	deque<Output_batch*> held;
	bool writer_busy;
	bool writer_stopping;
	bool flush_requested;

	//The batch check of the other side's norms (on the GPU): a Batch_smooth of the other side, with the primes up to OTHER_STAGE_BOUND
	//and two large primes, which the writer makes when the first norms come. A relation was already found, so its norm is smooth: the
	//tree leaves a cofactor of the large primes, which the checks split, and the rest is factored with the leaves. Its results come
	//back to the writer by number (ticket), through result_sink. For the norms with more large primes (about a quarter on SNFS-240),
	//other_factoring factors the cofactor and the rest of the norm apart (or the whole norm, for the few with a cofactor of more than
	//OTHER_STAGE_MAX_COFACTOR_BITS bits), in jobs for each tree of other_stage, which the writer waits for when it writes them. Most of
	//the primes of the norms are then found by the tree instead of Pollard rho.
	Batch_smooth* other_stage;
	bool other_stage_tried, other_stage_finished;
	long long other_tickets; //norms given to other_stage
	enum other_state { OTHER_FACTORED, OTHER_FAILED, OTHER_IN_JOB };
	struct Other_result
	{
		other_state state;
		vector<unsigned long long> primes; //all prime factors of the norm, in increasing order (OTHER_FACTORED)
		int job; //of other_factoring (OTHER_IN_JOB)
	};
	unordered_map<long long, Other_result> other_results;
	vector<long long> fallback_tickets, fallback_a, fallback_b; //not factored by other_stage, for the next job
	vector<long long> split_tickets; //their cofactor and the rest factored apart, for the next job of values:
	vector<unsigned long long> split_limbs; //the rest and the cofactor of each, at FACTOR_LIMBS limbs each
	vector<unsigned int> split_sizes;
	unordered_map<int, vector<long long>> job_tickets; //the tickets of the jobs not yet taken, in their order in the job
	unordered_map<int, bool> job_split; //a job of values, two for each ticket
	deque<int> pending_jobs; //in the order they were submitted, which is the order their results must be taken in
	void Create_other_stage();
	void Submit_fallback();
	void Take_job(int job);
	void Take_one_job(int job);
	//Set if this is an other_stage. For a norm that was not factored, cofactor is the part with the primes above the bound of the tree,
	//if it is a candidate (NULL if not), and norm the norm.
	function<void(long long ticket, long long a, long long b, bool factored, vector<unsigned long long>& primes, const mpz_class* cofactor, const mpz_class* norm)> result_sink;
	long long sink_next; //the ticket of the next relation of an other_stage
	mutex output_lock;
	condition_variable output_changed;
	thread writer_thread;

	string save_filename;
	FILE* output; //open from the first batch on; flushed after each batch, so that its size is that of the relations written
	int tree_size;
	int max_entries_before_check;
	atomic<long long> num_relations_found; //when the relations of a batch are checked, less those whose other side fails later
	atomic<long long> num_relations_written; //by this object as a writer
	atomic<long long> num_relations_checked;
	atomic<long long> num_factoring_failures;
	atomic<long long> num_other_side_failures;
	int num_relations_per_batch;
	int min_saved_prime_limit;

	Polynomial* poly;
	Polynomial* other_poly; //for rebuilding the other side's factor list
	bool rebuild_all; //for every relation, not only those without a saved list
public:
	//"format: keys": write only a,b of the relations found, without factoring them or rebuilding the other side (for building
	//a store of pairs, where only smoothness matters).
	bool keys_only = false;
private:

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
	function<void(const long long position[2], long long output_size, long long written)> batch_done;
	long long output_size; //of the output file after the last batch written
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
