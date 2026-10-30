#pragma once
#include "mpirxx.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include "Polynomial.h"
#include "Fft_product.h"
using namespace std;

#define MAX_SAVED_PRIMES 64 //Most primes in the saved factor list of a relation from a binary relation file.
#define MAX_COFACTOR_PRIMES 96 //Most prime factors of a 3LP cofactor, which is less than 2^96.

class Batch_smooth
{
public:
	//other_side_poly: if not NULL, the factor list of the side that is not batch checked is rebuilt from its norm, with all primes and their multiplicities (as CADO-NFS requires).
	//If NULL, that factor list is printed as it was read.
	//three_lp_stage: if true, this is the 3LP second stage. Its inputs are the cofactors given with Add_candidate, which are checked against all primes below 2^lbp (fblim and mfb are not used).
	Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side, int print_limit, Polynomial* other_side_poly, int max_tree_levels, bool three_lp_stage = false);
	~Batch_smooth(void);
	void Do_batch_check();
	void Finish(); //Check the remaining relations, and any 3LP candidates in the second stage. Call once, after the last relation.
	bool Add_number(long long a_val, long long b_val, const char* relation_factor_list);
	bool Add_number_binary(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes);
	bool Add_candidate(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes, mpz_srcptr cofactor);
	long long Get_num_relations_found();

	//3LP candidates, cofactors from 2 lbp + 1 to mfb bits, are passed to a second stage that this object then owns (NULL: not used),
	//and/or appended to a candidate file. Returns false if the file could not be opened.
	bool Set_three_lp_handling(Batch_smooth* second_stage, string candidate_filename);

	int side;

private:
	mpz_class* prime_product;
	mpz_class** inputs;
	long long * a_vals, * b_vals;
	bool* is_smooth; // set by the batch check for relations that are smooth.
	char** list_of_factors; // holds the already known smooth side factor list, as text or (binary_lists) as the number of primes followed by the primes.
	bool binary_lists;
	void Get_saved_list(int index, string& list);
	void Get_saved_primes(int index, vector<unsigned long long>& primes);
	mpz_class* temp_number;
	mpz_class* relation_value;

	int num_factors;
	unsigned int* relation_factor_list; //for Pollard rho factoring of the found relations.
	unsigned int* other_side_factor_list; //for factoring the norm of the other side, when rebuilding its factor list.
	unsigned long long* lp1;
	unsigned long long* lp2;
	unsigned int large_primes[MAX_COFACTOR_PRIMES]; //the large primes of the relation being printed.
	int num_large_primes;

	void Setup_prime_product(long long smooth_bound);
	void Setup_memory(int tree_levels, int max_bits_in_input);

	void Output_solutions();
	bool Complete_other_side(int index, string& list);
	bool Timed_complete_other_side(int index, string& list);

	string save_filename;
	int tree_size;
	int entries;
	int start_entry;
	int max_entries_before_check;
	long long num_relations_found;
	long long num_factoring_failures;
	long long num_other_side_failures;
	int num_relations_per_batch;
	int min_saved_prime_limit;

	Polynomial* poly;
	Polynomial* other_poly; // NULL unless the other side's factor list should be rebuilt.

	bool using_snfs_deg4_binomial;
	bool binomial_fits; //the coefficients of the binomial fit in 63 bits, for Binomial_value_native.
	long long binomial_c4, binomial_c0;
	bool Binomial_value_native(mpz_ptr result, long long a, long long b);
	long long max_batch_prime;

	int one_lp_bit_limit;
	long long two_lp_min_limit;
	int max_cofactor_bit_size;
	bool do_cofactorization;
	bool check_two_lp_smooth(mpz_class& number, unsigned int index);
	mpz_class *factor1, *factor2;

	//3LP handling.
	bool three_lp_stage; //this object is the second stage.
	int three_lp_max_bits; //largest 3LP candidate, mfb if it is more than 2 lbp, otherwise 0 (no 3LP).
	mpz_class three_lp_min; //smallest 3LP candidate, fblim^3.
	vector<int> candidate_index; //3LP candidates of the current batch, sent on after the batch is checked.
	vector<mpz_class> candidate_cofactor;
	Batch_smooth* second_stage;
	FILE* candidate_file;
	long long num_candidates;
	long long num_candidates_skipped; //could not be saved (a, b or saved primes too large for the candidate format).
	void Run_batch_check();
	void Send_candidates();

	void Scaled_root_value(mpz_ptr root, mp_bitcnt_t precision);
	void Barrett_reduce(mpz_ptr r, mpz_srcptr a, mpz_srcptr x, mp_bitcnt_t e, mp_bitcnt_t bits);
	mpz_class root_shifted, root_inverse, root_remainder, root_quotient, root_temp, root_part; //temporaries for Scaled_root_value.
	Cyclic_product cyclic; //wraparound products of the remainder tree.
	mpz_class wrap_result[2];
};

