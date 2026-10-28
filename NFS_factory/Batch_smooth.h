#pragma once
#include "mpirxx.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include "Polynomial.h"
using namespace std;

#define MAX_SAVED_PRIMES 64 //Most primes in the saved factor list of a relation from a binary relation file.

class Batch_smooth
{
public:
	//other_side_poly: if not NULL, the factor list of the side that is not batch checked is rebuilt from its norm, with all primes and their multiplicities (as CADO-NFS requires).
	//If NULL, that factor list is printed as it was read.
	Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side, int print_limit, Polynomial* other_side_poly, int max_tree_levels);
	~Batch_smooth(void);
	void Do_batch_check();
	bool Add_number(long long a_val, long long b_val, const char* relation_factor_list);
	bool Add_number_binary(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes);
	int Get_num_relations_found();

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
	long num_relations_found;
	long long num_factoring_failures;
	long long num_other_side_failures;
	int num_relations_per_batch;
	int min_saved_prime_limit;

	Polynomial* poly;
	Polynomial* other_poly; // NULL unless the other side's factor list should be rebuilt.

	bool using_snfs_deg4_binomial;
	long long max_batch_prime;

	int one_lp_bit_limit;
	long long two_lp_min_limit;
	int max_cofactor_bit_size;
	bool do_cofactorization;
	bool check_two_lp_smooth(mpz_class& number, unsigned int index);
	mpz_class *factor1, *factor2;

	void Middle_product(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t shift, mp_bitcnt_t precision);
	void Scaled_root_value(mpz_ptr root, mp_bitcnt_t precision);
	mpz_class root_shifted, root_inverse, root_remainder, root_quotient; //temporaries for Scaled_root_value.
	vector<mp_limb_t> wrap_buffer; //scratch space for the wraparound products.
	mpz_class wrap_result[2];
};

