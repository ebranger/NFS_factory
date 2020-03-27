#pragma once
#include "mpirxx.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include "Polynomial.h"
using namespace std;

class Batch_smooth
{
public:
	Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side, int print_limit);
	~Batch_smooth(void);
	void Do_batch_check();
	void Add_number( long long a_val, long long b_val, char* relation_factor_list);
	int Get_num_relations_found();

	int side;

private:
	mpz_class* prime_product;
	mpz_class** inputs;
	long long * a_vals, * b_vals;
	char** list_of_factors; // holds the already known smooth side factor list.
	mpz_class* temp_number;
	mpz_class* relation_value;
	mpz_class* small_factors;

	int num_factors;
	unsigned int* relation_factor_list; //for Pollard rho factoring of the found relations.
	unsigned long long* lp1;
	unsigned long long* lp2;

	bool product_set_up;
	void Setup_prime_product(int smooth_bound);
	void Setup_memory(int tree_levels, int max_bits_in_input);
	
	void Output_solutions();

	string save_filename;
	int tree_size;
	int entries;
	int start_entry;
	int max_entries_before_check;
	long num_relations_found;
	int num_relations_per_batch;
	int min_saved_prime_limit;

	Polynomial* poly;

	bool using_snfs_deg4_binomial;
	int max_batch_prime;
	
	int one_lp_bit_limit;
	long long two_lp_min_limit;
	int max_cofactor_bit_size;
	bool do_cofactorization;
	bool check_two_lp_smooth(mpz_class number, unsigned int index);
	int smooth_factor_reduce;
	mpz_class *factor1, *factor2;
};

