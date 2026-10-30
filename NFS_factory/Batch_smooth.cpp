#include "Batch_smooth.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <cerrno>
#include <time.h>
#include <limits.h>
#include "cofactorize.h"

#include "Profile.h"
#include "Prime_product.h"
#include "Smooth_factor.h"

using namespace std;

#define MAX_TREE_LEVELS 20 //Default for the maximum number of levels, can be changed with treelevels in the polynomial file.
//max 22 requires ~3GB of RAM. However, for SNFS-220 I seem to get the same speed for lower size, which saves a lot of RAM. 20 is ~7% slower, but saved 2GB of RAM, so worth it.
//For larger factorizations, consider increasing this a little, but for every +1 the memory requirement for the tree roughly doubles.

#define MAX_SQUARINGS_IN_CHECK 1 //Bernsteins algorithm does squaring at the bottom nodes, to handle powers of primes. This is the largest power checked.
//may want to do a small-prime check spearately, so that only the larger primes are still present when checking.
// Note: max squarings 2 means max power is 2^2=4.
//Currently: do only 1 squaring, it misses a few very rare relations but increase speed slightly.

#define MAX_CHARS_IN_RELATION_LIST 320 //Longest saved factor list. Relations with longer lists are rejected. Typical lists are less than 100 characters.

#define MAX_FACTORS_IN_RELATION 512 //Each factor is at least 2, so this is enough for any relation value of up to 512 bits.

#define WRAPAROUND_MIN_LIMBS 512 //Use wraparound products in the remainder tree for nodes of at least this many limbs. Tuned on an i7-14700F.

static void append_prime(string& list, unsigned long long prime);

//Initialize the batch smoothness class
Batch_smooth::Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side, int print_limit, Polynomial* other_side_poly, int max_tree_levels, bool three_lp_stage)
{
	//init.
	temp_number = new mpz_class;

	long long max_prime_limit = 1LL << lbp;

	this->three_lp_stage = three_lp_stage;
	if (three_lp_stage)
	{
		//The cofactors have no prime factors up to the fblim of the first stage, but the primes below that are cheap to include,
		//and then the candidates do not depend on the fblim used to find them.
		fblim = max_prime_limit - 1;
		mfb = 3 * lbp;
	}

	//if fblim and 2^lbp match, then we batch-check all numbers against the full primorial and do no cofactorization of any remaining numbers.
	// fblim > 2^lbp makes little sense? so assume equal, otherwise we will end up printing relations with some factors too large. Maybe check this in input file and warn user?
	if (max_prime_limit >= fblim && !three_lp_stage)
	{
		do_cofactorization = true;
	}
	else
	{ 
		do_cofactorization = false;
	}

	poly = input_poly;
	other_poly = other_side_poly;
	binary_lists = false;

	if (poly->poly_is_binomial_SNFS_deg4())
	{
		using_snfs_deg4_binomial = true;
	}
	else
	{
		using_snfs_deg4_binomial = false;
	}

	//The coefficients of a binomial, for the native evaluation, if they fit in 63 bits.
	binomial_fits = using_snfs_deg4_binomial && mpz_cmpabs_ui(poly->get_coefficient(4).get_mpz_t(), (1ULL << 63) - 1) <= 0 && mpz_cmpabs_ui(poly->get_coefficient(0).get_mpz_t(), (1ULL << 63) - 1) <= 0;
	binomial_c4 = binomial_fits ? mpz_get_si(poly->get_coefficient(4).get_mpz_t()) : 0;
	binomial_c0 = binomial_fits ? mpz_get_si(poly->get_coefficient(0).get_mpz_t()) : 0;

	side = batch_checking_side;

	prime_product = new mpz_class();
	*prime_product = 1;

	Setup_prime_product(fblim);
	max_batch_prime = fblim;

	//For handling 2lp cofactorization.
	one_lp_bit_limit = lbp;
	two_lp_min_limit = 2 * (long long)floor(log2((double)fblim));
	max_cofactor_bit_size = min(mfb, 2 * lbp);

	//Larger cofactors, up to mfb bits, can have three large primes. They must be composite and at least fblim^3.
	three_lp_max_bits = (mfb > 2 * lbp && !three_lp_stage) ? mfb : 0;
	mpz_ui_pow_ui(three_lp_min.get_mpz_t(), (mpir_ui)fblim, 3);
	second_stage = NULL;
	candidate_file = NULL;
	num_candidates = 0;
	num_candidates_skipped = 0;

	entries = 0;

	//Estiamte a good tree size for the remainder tree
	//Choose the smallest tree where the product of the relation values at the root is at least half the size of the prime product.
	//A larger tree makes the prime product modulo the root cheaper, but costs one more level of products and twice the memory.
	//For SNFS-240 on an i7-14700F, the root being 0.57 times the prime product (19 levels) was as fast as 1.15 times (20 levels),
	//and 0.29 times (18 levels) about 5 % slower.
	tree_size = (int)ceil(log2((double)mpz_sizeinbase(prime_product->get_mpz_t(), 2) / max_bits_in_input));
	tree_size = min((max_tree_levels > 0) ? max_tree_levels : MAX_TREE_LEVELS, tree_size);  // Force the remainder tree to not be too big to save memory.
	tree_size = max(2, tree_size); // The remainder tree needs at least one level above the relations.

	start_entry = (1 << (tree_size - 1)) - 1;
	max_entries_before_check = 1 << (tree_size - 1);
	save_filename = filename;
	num_relations_found = 0;
	num_factoring_failures = 0;
	num_other_side_failures = 0;
	min_saved_prime_limit = print_limit; //Do not print primes less than this. Msieve reconstructs primes less than 1000, CADO needs all primes (so set this to 1 to print everything).

	Setup_memory(tree_size, max_bits_in_input);

	relation_factor_list = new unsigned int[MAX_FACTORS_IN_RELATION]; //Small primes are stored with multiplicity, so this needs room for many factors.
	other_side_factor_list = new unsigned int[MAX_FACTORS_IN_RELATION];

	factor1 = new mpz_class;
	factor2 = new mpz_class;
}

//Destructor, free all memory that was requested. Relations must already have been checked with Do_batch_check.
Batch_smooth::~Batch_smooth(void)
{
	if (num_factoring_failures > 0)
	{
		cout << "Could not fully factor " << num_factoring_failures << " smooth relations, they were not saved." << endl;
	}

	if (num_other_side_failures > 0)
	{
		cout << num_other_side_failures << " smooth relations had a factor list for the other side that did not match its norm, they were not saved." << endl;
	}

	delete second_stage;
	if (candidate_file != NULL)
	{
		fclose(candidate_file);
	}

	for (int i = 0; i < 2 * num_relations_per_batch - 1; i++)
	{
		delete inputs[i];
	}
	delete[] inputs;
	delete prime_product;
	delete temp_number;

	for (int i = 0; i < num_relations_per_batch; i++)
	{
		delete[] list_of_factors[i];
	}
	delete[] list_of_factors;
	delete[] a_vals;
	delete[] b_vals;
	delete[] is_smooth;
	delete[] lp1;
	delete[] lp2;
	delete[] relation_factor_list;
	delete[] other_side_factor_list;
	delete factor1;
	delete factor2;
}

//Setup the primorials to be used.
void Batch_smooth::Setup_prime_product(long long smooth_bound)
{
	//So apparently mpir knows primorials. Problem solved!

	if (smooth_bound < 1000)
	{
		//I expect a minimum bound of ~2^20, anything less is probably a bad idea. 
		//Use 1000 as a basic sanity check. 
		return;
	}

	//Only the primes that can divide the norms are included, with small primes to higher powers.
	setup_prime_product(poly, (uint32)smooth_bound, prime_product);

	size_t size = mpz_sizeinbase(prime_product->get_mpz_t(), 2);
	cout << "Size of prime product: " << size << " bits." << endl;

	return;
}

//allocate sufficient memory to hold the remainder trees. Since their size may change from run to run, allocate enough space that no realloc needs to be done. 
//to check: is the MPIR mpz_class memory mangement efficient enough that this can be skipped?
void Batch_smooth::Setup_memory(int tree_levels, int max_bits_in_input)
{
	//allocate enough memory to hold all the numbers in the product tree without having to increase any mpz size
	int temp = 1 << tree_levels;
	inputs = new mpz_class*[temp - 1];
	lp1 = new unsigned long long[temp];
	lp2 = new unsigned long long[temp];
	int limit;

	temp = temp / 2;
	num_relations_per_batch = temp;
	a_vals = new long long[num_relations_per_batch];
	b_vals = new long long[num_relations_per_batch];
	is_smooth = new bool[num_relations_per_batch];

	list_of_factors = new char*[num_relations_per_batch];

	for (int i = 0; i < temp; i++)
	{
		//allocate space for storing the relation factor list
		list_of_factors[i] = new char[MAX_CHARS_IN_RELATION_LIST];

		//set large prime list to 0
		lp1[i] = 0;
		lp2[i] = 0;
		is_smooth[i] = false;
	}

	temp = max_bits_in_input;

	for (int level = tree_levels; level >= 1; level--)
	{
		limit = (1 << (level - 1)) - 1;

		for (int i = limit; i <= 2 * limit; i++)
		{
			//allocate an appropriate amount of memory, to aviod having to grow the mpz.
			inputs[i] = new mpz_class();
			mpz_init2(inputs[i]->get_mpz_t(), temp);
			*inputs[i] = 1; //multiplicative unity.
		}
		temp = temp * 2;
	}
}

//Add one relation to be batch-smoothness checked. Returns false if the factor list is too long to be saved.
bool Batch_smooth::Add_number(long long a_val, long long b_val, const char* relation_factor_list)
{
	//add one relation for batch-factoring.
	if (strlen(relation_factor_list) >= MAX_CHARS_IN_RELATION_LIST)
	{
		return false;
	}

	a_vals[entries] = a_val;
	b_vals[entries] = b_val;
	is_smooth[entries] = true; //Until the batch check shows otherwise.
	strcpy(list_of_factors[entries], relation_factor_list);
	binary_lists = false;

	entries++;

	//Is the tree full?
	if (entries == max_entries_before_check) //entries points at first empty position.
	{
		Do_batch_check();
	}
	return true;
}

//Add one relation to be batch-smoothness checked, with the factor list of the other side given as primes, as read from a binary relation file.
//The primes are only converted to text if the relation is smooth. Returns false if there are too many primes to be saved.
bool Batch_smooth::Add_number_binary(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes)
{
	if (num_primes > MAX_SAVED_PRIMES)
	{
		return false;
	}

	a_vals[entries] = a_val;
	b_vals[entries] = b_val;
	is_smooth[entries] = true; //Until the batch check shows otherwise.

	//The saved list holds the number of primes, followed by the primes.
	unsigned int* saved = reinterpret_cast<unsigned int*>(list_of_factors[entries]);
	saved[0] = num_primes;
	memcpy(saved + 1, primes, num_primes * sizeof(unsigned int));
	binary_lists = true;

	entries++;

	//Is the tree full?
	if (entries == max_entries_before_check)
	{
		Do_batch_check();
	}
	return true;
}

//Add one 3LP candidate to the second stage: a relation with the factor list of the other side, and the cofactor of the batch checked side
//that is left after dividing out the primes up to fblim of the first stage. Returns false if there are too many primes to be saved.
bool Batch_smooth::Add_candidate(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes, mpz_srcptr cofactor)
{
	if (num_primes > MAX_SAVED_PRIMES)
	{
		return false;
	}

	//The cofactor is checked instead of the polynomial value.
	mpz_set(inputs[start_entry + entries]->get_mpz_t(), cofactor);
	return Add_number_binary(a_val, b_val, primes, num_primes);
}

bool Batch_smooth::Set_three_lp_handling(Batch_smooth* second_stage, string candidate_filename)
{
	this->second_stage = second_stage;

	if (!candidate_filename.empty())
	{
		candidate_file = fopen(candidate_filename.c_str(), "ab");
		if (candidate_file == NULL)
		{
			return false;
		}
	}
	return true;
}

//Pass on the 3LP candidates of the batch that was just checked, to the second stage and/or the candidate file.
//The candidate file has the binary relation format with the cofactor added: (int)a, (uint)b, (uint)num_primes, (uint)prime1 ... (uint)prime_num_primes,
//(uint64)low and (uint64)high 64 bits of the cofactor.
void Batch_smooth::Send_candidates()
{
	vector<unsigned long long> primes;
	unsigned int saved[MAX_SAVED_PRIMES];

	for (size_t n = 0; n < candidate_index.size(); n++)
	{
		int index = candidate_index[n];
		Get_saved_primes(index, primes);

		bool fits = primes.size() <= MAX_SAVED_PRIMES && a_vals[index] >= INT_MIN && a_vals[index] <= INT_MAX && b_vals[index] <= UINT_MAX;
		for (size_t j = 0; j < primes.size() && fits; j++)
		{
			fits = primes[j] <= UINT_MAX;
			saved[j] = (unsigned int)primes[j];
		}

		if (!fits)
		{
			num_candidates_skipped++;
			continue;
		}

		num_candidates++;
		unsigned int num_primes = (unsigned int)primes.size();

		if (candidate_file != NULL)
		{
			int a = (int)a_vals[index];
			unsigned int b = (unsigned int)b_vals[index];
			unsigned long long cofactor[2] = { mpz_getlimbn(candidate_cofactor[n].get_mpz_t(), 0), mpz_getlimbn(candidate_cofactor[n].get_mpz_t(), 1) };
			fwrite(&a, sizeof(int), 1, candidate_file);
			fwrite(&b, sizeof(unsigned int), 1, candidate_file);
			fwrite(&num_primes, sizeof(unsigned int), 1, candidate_file);
			fwrite(saved, sizeof(unsigned int), num_primes, candidate_file);
			fwrite(cofactor, sizeof(unsigned long long), 2, candidate_file);
		}

		if (second_stage != NULL)
		{
			second_stage->Add_candidate(a_vals[index], b_vals[index], saved, num_primes, candidate_cofactor[n].get_mpz_t());
		}
	}

	candidate_index.clear();
}

void Batch_smooth::Finish()
{
	Do_batch_check();

	if (three_lp_max_bits > 0)
	{
		cout << "Found " << num_candidates << " 3LP candidates";
		if (num_candidates_skipped > 0)
		{
			cout << " (and skipped " << num_candidates_skipped << " with a, b or saved primes too large for the candidate format)";
		}
		cout << "." << endl;
	}

	if (second_stage != NULL)
	{
		second_stage->Finish();
		cout << "The 3LP second stage found " << second_stage->Get_num_relations_found() << " relations." << endl;
	}
}

//The saved factor list of the other side of a relation, as text.
void Batch_smooth::Get_saved_list(int index, string& list)
{
	if (!binary_lists)
	{
		list = list_of_factors[index];
		return;
	}

	const unsigned int* saved = reinterpret_cast<const unsigned int*>(list_of_factors[index]);
	list.clear();
	for (unsigned int i = 1; i <= saved[0]; i++)
	{
		append_prime(list, saved[i]);
	}
}

//The saved factor list of the other side of a relation, as numbers.
void Batch_smooth::Get_saved_primes(int index, vector<unsigned long long>& primes)
{
	primes.clear();

	if (binary_lists)
	{
		const unsigned int* saved = reinterpret_cast<const unsigned int*>(list_of_factors[index]);
		for (unsigned int i = 1; i <= saved[0]; i++)
		{
			primes.push_back(saved[i]);
		}
		return;
	}

	const char* text = list_of_factors[index];
	while (*text != '\0')
	{
		char* end;
		primes.push_back(strtoull(text, &end, 16));
		text = (*end == ',') ? end + 1 : end;
	}
}

//Check whether a remainder factors into two numbers less than the lp limit.
//If it does, it is a valid relation!
bool Batch_smooth::check_two_lp_smooth(mpz_class& number, unsigned int index)
{
	int cofactor_bits = (int)mpz_sizeinbase(number.get_mpz_t(), 2);

	if (cofactor_bits > max_cofactor_bit_size)
	{
		//Too big to be split into two primes both less than lp_limit. It can still have three, if it is composite and at least fblim^3.
		if (cofactor_bits <= three_lp_max_bits && number >= three_lp_min && !is_probable_prime(number.get_mpz_t()))
		{
			if (candidate_cofactor.size() <= candidate_index.size())
			{
				candidate_cofactor.emplace_back();
			}
			candidate_cofactor[candidate_index.size()] = number;
			candidate_index.push_back(index);
		}
		return false;
	}

	if (cofactor_bits > two_lp_min_limit)
	{
		//Can potentially split into two large primes, small enough to be useful
		if (is_probable_prime(number.get_mpz_t()))
		{
			return false;
		}

		*factor1 = 0;
		*factor2 = 0;

		if (cofactor_bits <= 63 && mpz_odd_p(number.get_mpz_t()))
		{
			//Pollard rho with 64-bit arithmetic.
			unsigned long long factor;
			if (split_composite64(mpz_get_ui(number.get_mpz_t()), &factor))
			{
				mpz_set_ui(factor1->get_mpz_t(), factor);
				mpz_divexact_ui(factor2->get_mpz_t(), number.get_mpz_t(), factor);
			}
		}

		if (*factor1 == 0)
		{
			//either rho failed or the cofactor is too big for 64-bit rho, so do qs.
			int qs_sucess = tinyqs(number.get_mpz_t(), factor1->get_mpz_t(), factor2->get_mpz_t());

			if (qs_sucess == 0)
			{
				//unsuccesful
				return false;
			}

		}

		//check if the factors are too big.
		if (mpz_sizeinbase(factor1->get_mpz_t(), 2) > one_lp_bit_limit)
		{
			return false;
		}
		if (mpz_sizeinbase(factor2->get_mpz_t(), 2) > one_lp_bit_limit)
		{
			return false;
		}

		//check that the factors are actual large primes. If it turns out to be a small prime that occured to a high
		//power, it can cause problems later when it is assumed to be a prime larger than the rlim/alim.
		if (!is_probable_prime(factor1->get_mpz_t()))
		{
			return false;
		}
		if (!is_probable_prime(factor2->get_mpz_t()))
		{
			return false;
		}

		//save the two found factors, the smaller first. I do not expect three factors to ever happen, though that may cause problems here...
		lp1[index] = min(factor1->get_ui(), factor2->get_ui());
		lp2[index] = max(factor1->get_ui(), factor2->get_ui());

		return true;
	}

	if (cofactor_bits > one_lp_bit_limit)
	{
		//Should be a prime and too big.
		return false;
	}

	if (number == 1)
	{
		//no large primes found, finished!
		return true;
	}

	//should be one large prime by now. Make sure it is actually a large prime as well.
	if (number > max_batch_prime && is_probable_prime(number.get_mpz_t()))
	{
		lp1[index] = mpz_get_ui(number.get_mpz_t());
		lp2[index] = 0;
	}
	else
	{
		//I have seen this happen when the remainder was a product of one large prime and one very small that occured to a high power in the relation,
		//and the small prime was not found by the batch code due to the high power.
		//Since this seems to happend very rarely, for now just do a basic check if the sizes of the two factors makes the relation valid,
		//and let the relation factoring code handle it later, without any lp hints from here.
		num_factors = factor_general(relation_factor_list, number.get_mpz_t(), MAX_FACTORS_IN_RELATION);

		if (num_factors < 0)
		{
			return false;
		}

		for (int i = 0; i < num_factors; i++)
		{
			if ((unsigned long long)relation_factor_list[i] >> one_lp_bit_limit != 0)
			{
				//factor is 2^lbp or larger.
				return false;
			}
		}

	}

	//if the remainder after dividing out all small primes is less than the limit, it is prime and smooth (or possibly composite with high-power primes in relation, but still OK.)
	return true;
}

//|c| x^4 for |x| < 2^32 and |c| < 2^63, which is less than 2^191, as three limbs.
static void binomial_term(unsigned long long term[3], long long c, long long x)
{
	unsigned long long abs_x = (unsigned long long)(x < 0 ? -x : x);
	unsigned long long abs_c = (unsigned long long)(c < 0 ? -c : c);
	unsigned long long square = abs_x * abs_x;
	unsigned long long power_high;
	unsigned long long power_low = _umul128(square, square, &power_high);
	unsigned long long high0, high1;
	term[0] = _umul128(power_low, abs_c, &high0);
	unsigned long long middle = _umul128(power_high, abs_c, &high1);
	unsigned char carry = _addcarry_u64(0, middle, high0, &term[1]);
	_addcarry_u64(carry, high1, 0, &term[2]);
}

//The absolute value of the binomial c4 a^4 + c0 b^4 with native arithmetic. Returns false, and leaves the result alone, if the
//coefficients or a and b are too large for it.
bool Batch_smooth::Binomial_value_native(mpz_ptr result, long long a, long long b)
{
	const long long limit = 1LL << 32;
	if (!binomial_fits || a <= -limit || a >= limit || b <= -limit || b >= limit)
	{
		return false;
	}

	unsigned long long term4[3], term0[3], value[3];
	binomial_term(term4, binomial_c4, a);
	binomial_term(term0, binomial_c0, b);

	if ((binomial_c4 < 0) == (binomial_c0 < 0))
	{
		unsigned char carry = _addcarry_u64(0, term4[0], term0[0], &value[0]);
		carry = _addcarry_u64(carry, term4[1], term0[1], &value[1]);
		_addcarry_u64(carry, term4[2], term0[2], &value[2]);
	}
	else
	{
		//The difference of the terms, the smaller from the larger.
		unsigned long long* larger = term4;
		unsigned long long* smaller = term0;
		if (term0[2] > term4[2] || (term0[2] == term4[2] && (term0[1] > term4[1] || (term0[1] == term4[1] && term0[0] > term4[0]))))
		{
			larger = term0;
			smaller = term4;
		}
		unsigned char borrow = _subborrow_u64(0, larger[0], smaller[0], &value[0]);
		borrow = _subborrow_u64(borrow, larger[1], smaller[1], &value[1]);
		_subborrow_u64(borrow, larger[2], smaller[2], &value[2]);
	}

	if (result->_mp_alloc < 3)
	{
		mpz_realloc2(result, 3 * GMP_LIMB_BITS);
	}
	int size = 3;
	while (size > 0 && value[size - 1] == 0)
	{
		size--;
	}
	for (int i = 0; i < size; i++)
	{
		result->_mp_d[i] = value[i];
	}
	result->_mp_size = size;
	return true;
}

//Check the relations added so far. The time of the 3LP second stage is counted as one phase of the profile, with its phases kept separately.
void Batch_smooth::Do_batch_check()
{
	if (!three_lp_stage)
	{
		Run_batch_check();
		return;
	}

	unsigned long long outer_cycles[PROF_NUM_PHASES];
	memcpy(outer_cycles, profile_cycles, sizeof(outer_cycles));
	unsigned long long t = profile_now();

	Run_batch_check();

	for (int i = 0; i < PROF_NUM_PHASES; i++)
	{
		profile_cycles_3lp[i] += profile_cycles[i] - outer_cycles[i];
	}
	memcpy(profile_cycles, outer_cycles, sizeof(outer_cycles));
	profile_add(PROF_3LP, t);
}

//Perform the remainder tree calculations.
//Core parts of this function was taken from the CADO-NFS package.
//Use a scaled remainder tree for speed.
void Batch_smooth::Run_batch_check()
{
	int limit;

	if (entries == 0)
	{
		//Nothing to check.
		return;
	}

	unsigned long long t = profile_now();

	//set up polynomial values. The second stage checks the cofactors that were added.
	for (int index = 0; index < entries && !three_lp_stage; index++)
	{
		//use different poly evaluation if polynomial has form k*x^4+l which is faster.
		if (using_snfs_deg4_binomial)
		{
			if (!Binomial_value_native(inputs[start_entry + index]->get_mpz_t(), a_vals[index], b_vals[index]))
			{
				poly->evaluate_binomial_polynomial_homogeneous_SNFS_deg4(inputs[start_entry + index], a_vals[index], b_vals[index], temp_number);
			}
		}
		else
		{
			poly->evaluate_polynomial_homogenous(inputs[start_entry + index], a_vals[index], b_vals[index]);
		}

		if (*inputs[start_entry + index] < 0)
		{
			*inputs[start_entry + index] *= -1;
		}
	}

	t = profile_add(PROF_POLY_EVAL, t);

	//set up product tree.
	for (int level = tree_size; level > 1; level--)
	{
		limit = (1 << (level - 2)) - 1;

		for (int i = limit; i <= 2 * limit; i++)
		{
			*inputs[i] = *inputs[2 * i + 1] * *inputs[2 * i + 2];
		}
	}

	t = profile_add(PROF_PRODUCT_TREE, t);

	size_t size = mpz_sizeinbase(inputs[0]->get_mpz_t(), 2);
	cout << "Size of top node of product tree: " << size << " bits." << endl;


	//The core remainder tree calculations are based on the CADO-nfs batch function

	long h = tree_size + 2; // extra fudge bits.
	unsigned long i, guard, *product_bits;
	product_bits = new unsigned long[2 * num_relations_per_batch];

	guard = h;

	//Get the size of all products, they are needed for shifts later.
	for (i = 0; i < (unsigned long)(2 * num_relations_per_batch - 1); i++)
	{
		product_bits[i] = (unsigned long)mpz_sizeinbase(inputs[i]->get_mpz_t(), 2);
	}

	mpz_class Q;

	// now the root of the tree is prime_product / product_of_input_numbers * 2^(bits in product of input number + a few guard bits)
	Scaled_root_value(inputs[0]->get_mpz_t(), product_bits[0] + guard);
	t = profile_add(PROF_ROOT, t);

	for (int level = 2; level <= tree_size - 1; level++)
	{
		limit = (1 << (level - 1)) - 1;

		for (int i = limit; i <= 2 * limit; i = i + 2)
		{
			//calculate the next level of the tree.
			// so if the parent node is Num / (leaf1 * leaf2), multiply by the leaves and right-shift since the fraction is now smaller, requiring fewer bits of storage.
			unsigned long parent_bits = product_bits[(i - 1) / 2];

			if ((parent_bits + guard) / GMP_LIMB_BITS >= WRAPAROUND_MIN_LIMBS)
			{
				//Large numbers: only the middle part of the products is needed, which is faster to get with a wraparound product.
				cyclic.Middle_product_pair(wrap_result[0].get_mpz_t(), wrap_result[1].get_mpz_t(), inputs[(i - 1) / 2]->get_mpz_t(), inputs[i + 1]->get_mpz_t(), inputs[i]->get_mpz_t(),
					parent_bits - product_bits[i], parent_bits - product_bits[i + 1], parent_bits + guard);
				mpz_swap(inputs[i]->get_mpz_t(), wrap_result[0].get_mpz_t());
				mpz_swap(inputs[i + 1]->get_mpz_t(), wrap_result[1].get_mpz_t());
				continue;
			}

			*inputs[i] = *inputs[(i - 1) / 2] * *inputs[i];
			*inputs[i + 1] = *inputs[(i - 1) / 2] * *inputs[i + 1];
			mpz_swap(inputs[i]->get_mpz_t(), inputs[i + 1]->get_mpz_t());

			mpz_tdiv_r_2exp(inputs[i]->get_mpz_t(), inputs[i]->get_mpz_t(), product_bits[(i - 1) / 2] + guard);
			mpz_div_2exp(inputs[i]->get_mpz_t(), inputs[i]->get_mpz_t(), product_bits[(i - 1) / 2] - product_bits[i]);

			mpz_tdiv_r_2exp(inputs[i + 1]->get_mpz_t(), inputs[i + 1]->get_mpz_t(), product_bits[(i - 1) / 2] + guard);
			mpz_div_2exp(inputs[i + 1]->get_mpz_t(), inputs[i + 1]->get_mpz_t(), product_bits[(i - 1) / 2] - product_bits[i + 1]);
		}
	}

	t = profile_add(PROF_REMAINDER_TREE, t);
	unsigned long long cofactor_start = profile_cycles[PROF_COFACTOR];

	//for the bottom node, i.e. the relation numbers, store the Prime_product % number in Q, to keep the number (in input_number[]). Then it can be used later for squaring tests mod the number.
	for (i = num_relations_per_batch - 1; i < (unsigned long)(2 * num_relations_per_batch - 1); i = i + 2)
	{
		//One bottom node
		Q = *inputs[(i - 1) / 2] * *inputs[i + 1];

		mpz_tdiv_r_2exp(Q.get_mpz_t(), Q.get_mpz_t(), product_bits[(i - 1) / 2] + guard);
		mpz_div_2exp(Q.get_mpz_t(), Q.get_mpz_t(), product_bits[(i - 1) / 2] - product_bits[i]);

		Q = Q * *inputs[i];
		mpz_div_2exp(Q.get_mpz_t(), Q.get_mpz_t(), product_bits[i]);
		Q = Q + (1UL << guard) - 1UL;
		mpz_div_2exp(Q.get_mpz_t(), Q.get_mpz_t(), guard);

		//Q = Q % *inputs[i]; // should not be needed, but valid relations are lost if this is not here... Except when i test it again later, in which case it is magically not needed again...

		for (int j = 0; j < MAX_SQUARINGS_IN_CHECK; j++)
		{
			Q = (Q * Q) % *inputs[i];
		}

		if (do_cofactorization == true)
		{
			//divide out all factors from the tree from the relation and put the result in Q. 
			mpz_gcd(temp_number->get_mpz_t(), inputs[i]->get_mpz_t(), Q.get_mpz_t());
			mpz_divexact(Q.get_mpz_t(), inputs[i]->get_mpz_t(), temp_number->get_mpz_t());

			//Check if Q is 2lp smooth
			unsigned long long tc = profile_now();
			bool smooth = check_two_lp_smooth(Q, i - num_relations_per_batch + 1);
			profile_add(PROF_COFACTOR, tc);
			if (!smooth)
			{
				//not smooth.
				is_smooth[i - num_relations_per_batch + 1] = false;
			}
		}
		else
		{
			if (Q > 0)
			{
				//Do not do 2lp testing. If Q = 0 then the number is smooth

				//not smooth.
				is_smooth[i - num_relations_per_batch + 1] = false;
			}
		}

		//The same for the other bottom node with the same parent node.
		Q = *inputs[(i - 1) / 2] * *inputs[i];

		mpz_tdiv_r_2exp(Q.get_mpz_t(), Q.get_mpz_t(), product_bits[(i - 1) / 2] + guard);
		mpz_div_2exp(Q.get_mpz_t(), Q.get_mpz_t(), product_bits[(i - 1) / 2] - product_bits[i + 1]);

		Q = Q * *inputs[i + 1];
		mpz_div_2exp(Q.get_mpz_t(), Q.get_mpz_t(), product_bits[i + 1]);
		Q = Q + (1UL << guard) - 1UL;
		mpz_div_2exp(Q.get_mpz_t(), Q.get_mpz_t(), guard);

		//Q = Q % *inputs[i+1]; // should not be needed, but valid relations are lost if this is not here... Except when i test it again later, in which case it is magically not needed again...


		for (int j = 0; j < MAX_SQUARINGS_IN_CHECK; j++)
		{
			Q = (Q * Q) % *inputs[i + 1];
		}

		if (do_cofactorization == true)
		{
			//divide out all factors from the tree from the relation and put the result in Q. 
			mpz_gcd(temp_number->get_mpz_t(), inputs[i + 1]->get_mpz_t(), Q.get_mpz_t());
			mpz_divexact(Q.get_mpz_t(), inputs[i + 1]->get_mpz_t(), temp_number->get_mpz_t());

			//Check if Q is 2lp smooth
			unsigned long long tc = profile_now();
			bool smooth = check_two_lp_smooth(Q, i - num_relations_per_batch + 2);
			profile_add(PROF_COFACTOR, tc);
			if (!smooth)
			{
				//not smooth.
				is_smooth[i - num_relations_per_batch + 2] = false;
			}
		}
		else
		{
			if (Q > 0)
			{
				//Do not do 2lp testing. If Q = 0 then the relation is smooth

				//not smooth.
				is_smooth[i - num_relations_per_batch + 2] = false;
			}
		}
	}

	profile_add(PROF_LEAVES, t);
	profile_cycles[PROF_LEAVES] -= profile_cycles[PROF_COFACTOR] - cofactor_start; //cofactorization is counted separately.

	Send_candidates();
	Output_solutions();

	delete[] product_bits;
	return;
}

//result = the bits of a from position up to position + count, as a number.
static void bits_of(mpz_ptr result, mpz_srcptr a, mp_bitcnt_t position, mp_bitcnt_t count)
{
	mp_size_t a_limbs = abs(a->_mp_size);
	mp_size_t first = (mp_size_t)(position / GMP_LIMB_BITS);
	if (count == 0 || first >= a_limbs)
	{
		mpz_set_ui(result, 0);
		return;
	}

	mp_size_t out_limbs = (mp_size_t)((count + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS);
	mp_size_t source_limbs = min(out_limbs + 1, a_limbs - first);
	if (result->_mp_alloc < source_limbs)
	{
		mpz_realloc2(result, source_limbs * GMP_LIMB_BITS);
	}

	if (position % GMP_LIMB_BITS != 0)
	{
		mpn_rshift(result->_mp_d, a->_mp_d + first, source_limbs, (unsigned int)(position % GMP_LIMB_BITS));
	}
	else
	{
		mpn_copyi(result->_mp_d, a->_mp_d + first, source_limbs);
	}

	out_limbs = min(out_limbs, source_limbs);
	if (count % GMP_LIMB_BITS != 0 && out_limbs == (mp_size_t)((count + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS))
	{
		result->_mp_d[out_limbs - 1] &= (((mp_limb_t)1) << (count % GMP_LIMB_BITS)) - 1;
	}
	while (out_limbs > 0 && result->_mp_d[out_limbs - 1] == 0)
	{
		out_limbs--;
	}
	result->_mp_size = (int)out_limbs;
}

//result = a modulo 2^L - 1, from 0 to 2^L - 2. temp is used as scratch space.
static void mod_2expm1(mpz_ptr result, mpz_srcptr a, mp_bitcnt_t L, mpz_ptr temp)
{
	mpz_set(result, a);
	while (mpz_sizeinbase(result, 2) > L)
	{
		mpz_tdiv_q_2exp(temp, result, L);
		mpz_tdiv_r_2exp(result, result, L);
		mpz_add(result, result, temp);
	}
	mpz_add_ui(result, result, 1);
	if (mpz_sizeinbase(result, 2) > L)
	{
		mpz_set_ui(result, 0);
	}
	else
	{
		mpz_sub_ui(result, result, 1);
	}
}

//Replace the root of the product tree X by the root of the scaled remainder tree, floor((P mod X) * 2^precision / X),
//where P is the prime product. One inverse of X is computed and used both for P mod X (Barrett reduction) and for the
//scaling, which is faster than two divisions that each compute an inverse. The result can be one too small, which the
//guard bits absorb, just as they absorb the truncation errors further down in the tree.
void Batch_smooth::Scaled_root_value(mpz_ptr root, mp_bitcnt_t precision)
{
	mp_bitcnt_t bits = mpz_sizeinbase(root, 2);
	mp_size_t limbs = abs(root->_mp_size);
	mp_bitcnt_t normalize = limbs * GMP_LIMB_BITS - bits;

	//V is the inverse of the root shifted to fill limbs + 1 limbs, root * V < 2^e, and at most INVERSE_MAX_ERROR less than
	//2^e / root, with e = 128 * limbs + 64 - normalize. The extra limb makes the scaling below accurate to less than one unit.
	mpz_mul_2exp(root_shifted.get_mpz_t(), root, normalize + GMP_LIMB_BITS);
	mpz_ptr v = root_inverse.get_mpz_t();
	cyclic.Inverse(v, root_shifted.get_mpz_t());
	mp_bitcnt_t e = 128 * limbs + GMP_LIMB_BITS - normalize;

	//r = P mod X. Barrett reduction works for numbers less than 2^(2 * bits), which the choice of tree size gives. A larger prime
	//product, as in the second stage when its tree is not full, is reduced one part at a time from the top: the remainder so far,
	//followed by the next bits - 1 bits of P, is less than 2^(2 * bits - 1).
	mpz_ptr r = root_remainder.get_mpz_t();
	mpz_srcptr p = prime_product->get_mpz_t();
	mp_bitcnt_t product_size = mpz_sizeinbase(p, 2);

	if (product_size <= 2 * bits)
	{
		Barrett_reduce(r, p, root, e, bits);
	}
	else
	{
		mp_bitcnt_t position = product_size;
		mpz_set_ui(r, 0);
		while (position > 0)
		{
			mp_bitcnt_t part = (position > bits - 1) ? bits - 1 : position;
			position -= part;
			bits_of(root_part.get_mpz_t(), p, position, part);
			mpz_mul_2exp(r, r, part);
			mpz_add(r, r, root_part.get_mpz_t());
			Barrett_reduce(r, r, root, e, bits);
		}
	}

	//floor(r * 2^precision / X) = floor(r * V / 2^(e - precision)), or one less.
	mpz_mul(root, r, v);
	mpz_tdiv_q_2exp(root, root, e - precision);
}

//r = a mod x, for a less than 2^(2 * bits), with x of bits bits and its inverse V from Scaled_root_value (root_inverse), so that
//x * V < 2^e. r can be a.
void Batch_smooth::Barrett_reduce(mpz_ptr r, mpz_srcptr a, mpz_srcptr x, mp_bitcnt_t e, mp_bitcnt_t bits)
{
	mp_bitcnt_t a_size = mpz_sizeinbase(a, 2);
	if (a_size < bits)
	{
		mpz_set(r, a);
		return;
	}

	//mu = floor(2^(2 * bits) / x), or up to two less, without the bits that do not matter for the quotient. The quotient has
	//at most quotient_bits bits, and the bits of mu more than 64 below that change the estimate by less than one.
	//The quotient estimate is then at most 5 too small.
	mpz_ptr q = root_quotient.get_mpz_t();
	mpz_ptr t = root_temp.get_mpz_t();
	mp_bitcnt_t quotient_bits = a_size - bits + 1;
	mp_bitcnt_t cut = (bits > quotient_bits + GMP_LIMB_BITS) ? bits - quotient_bits - GMP_LIMB_BITS : 0;
	mpz_tdiv_q_2exp(q, root_inverse.get_mpz_t(), e - 2 * bits + cut);
	mpz_tdiv_q_2exp(t, a, bits - 1);
	mpz_mul(q, q, t);
	mpz_tdiv_q_2exp(q, q, bits + 1 - cut);

	//r = a - q x is less than 6 x, so it can be computed modulo 2^L - 1 > 6 x, which only needs q x modulo 2^L - 1.
	mp_bitcnt_t L = cyclic.Product_mod(q, q, x, bits + 3);
	mod_2expm1(t, a, L, root_part.get_mpz_t());
	mpz_sub(r, t, q);
	if (mpz_sgn(r) < 0)
	{
		mpz_set_ui(q, 0);
		mpz_setbit(q, L);
		mpz_sub_ui(q, q, 1);
		mpz_add(r, r, q);
	}
	while (mpz_cmp(r, x) >= 0)
	{
		mpz_sub(r, r, x);
	}
}

//Append a prime to a comma separated list of hexadecimal primes.
static void append_prime(string& list, unsigned long long prime)
{
	char hex[20];
	sprintf(hex, "%llX", prime);

	if (!list.empty())
	{
		list += ',';
	}
	list += hex;
}

//Rebuild the factor list of the side that was not batch checked, with all primes and their multiplicities, as CADO-NFS requires.
//The saved list is checked against the norm, and primes left out by the siever (typically small ones) are added.
//Returns false if a saved prime does not divide the norm, or if the norm could not be fully factored.
bool Batch_smooth::Complete_other_side(int index, string& list)
{
	mpz_class norm;
	other_poly->evaluate_polynomial_homogenous(&norm, a_vals[index], b_vals[index]);
	mpz_abs(norm.get_mpz_t(), norm.get_mpz_t());

	if (norm == 0)
	{
		return false;
	}

	vector<unsigned long long> primes;
	Get_saved_primes(index, primes);

	//Divide out the primes from the saved list, once per time they are listed.
	for (size_t j = 0; j < primes.size(); j++)
	{
		if (primes[j] < 2 || !mpz_divisible_ui_p(norm.get_mpz_t(), (mpir_ui)primes[j]))
		{
			return false;
		}
		mpz_divexact_ui(norm.get_mpz_t(), norm.get_mpz_t(), (mpir_ui)primes[j]);
	}

	//Factor whatever is left, which also finds any higher multiplicities of the listed primes.
	if (norm > 1)
	{
		int num_other_factors = factor_general(other_side_factor_list, norm.get_mpz_t(), MAX_FACTORS_IN_RELATION);

		if (num_other_factors < 0)
		{
			return false;
		}

		for (int j = 0; j < num_other_factors; j++)
		{
			primes.push_back(other_side_factor_list[j]);
		}
	}

	sort(primes.begin(), primes.end());

	list.clear();
	for (size_t j = 0; j < primes.size(); j++)
	{
		append_prime(list, primes[j]);
	}
	return true;
}

//Complete_other_side, with the time counted in the profile.
bool Batch_smooth::Timed_complete_other_side(int index, string& list)
{
	unsigned long long t = profile_now();
	bool result = Complete_other_side(index, list);
	profile_add(PROF_OTHER_SIDE, t);
	return result;
}

//For each relation passing the batch smoothnes checking, complete the relation and print it.
void Batch_smooth::Output_solutions()
{
	//Print all found relations to file
	int limit;
	limit = (1 << (tree_size - 1)) - 1;
	int index;

	string checked_side_list; //factors of the relation value that was batch checked.
	string other_side_list; //factors of the other side, if it is rebuilt.
	string relation_text;
	char ab_text[50];

	int num_good = 0;

	unsigned long long t_start = profile_now();
	unsigned long long factor_start = profile_cycles[PROF_FACTOR] + profile_cycles[PROF_OTHER_SIDE];

	FILE* outfile;

	outfile = fopen(save_filename.c_str(), "a");

	if (outfile == NULL)
	{
		cout << "Could not open output file for writing relations." << endl;
		return;
	}

	for (int i = limit; i <= 2 * limit; i++)
	{
		index = i - limit;

		//If the relation passed the smoothness test, and has been checked.
		if (index < entries && is_smooth[index])
		{
			relation_value = inputs[i];

			//negative values are handled by msieve, ignore here.
			if (*relation_value < 0)
			{
				*relation_value = *relation_value * -1;
			}

			//factor the relation. First, divide out the large primes, which are printed first.
			num_factors = 0;
			num_large_primes = 0;

			if (three_lp_stage)
			{
				//The value is the cofactor from the first stage, and all its primes are below 2^lbp. Those are the large primes,
				//and the rest of the polynomial value is factored as for other relations.
				unsigned long long tf = profile_now();
				num_large_primes = factor_smooth(large_primes, relation_value->get_mpz_t(), MAX_COFACTOR_PRIMES);
				poly->evaluate_polynomial_homogenous(temp_number, a_vals[index], b_vals[index]);
				mpz_abs(temp_number->get_mpz_t(), temp_number->get_mpz_t());

				if (num_large_primes < 0 || *temp_number == 0 || !mpz_divisible_p(temp_number->get_mpz_t(), relation_value->get_mpz_t()))
				{
					num_factors = -1;
				}
				else
				{
					mpz_divexact(relation_value->get_mpz_t(), temp_number->get_mpz_t(), relation_value->get_mpz_t());
				}
				profile_add(PROF_FACTOR, tf);
			}
			else
			{
				if (lp1[index] > 1)
				{
					large_primes[num_large_primes++] = (unsigned int)lp1[index];
					mpz_divexact_ui(relation_value->get_mpz_t(), relation_value->get_mpz_t(), lp1[index]);
				}

				if (lp2[index] > 1)
				{
					large_primes[num_large_primes++] = (unsigned int)lp2[index];
					mpz_divexact_ui(relation_value->get_mpz_t(), relation_value->get_mpz_t(), lp2[index]);
				}
			}

			//run Pollard rho factoring on the remaining part of the relation. Since it is smooth Pollard rho should be a decent choice.
			if (num_factors == 0 && *relation_value > 1)
			{
				unsigned long long tf = profile_now();
				num_factors = factor_general(relation_factor_list, relation_value->get_mpz_t(), MAX_FACTORS_IN_RELATION);
				profile_add(PROF_FACTOR, tf);
			}

			if (num_factors < 0)
			{
				//Typically: one factor is too big for an uint32. An incomplete relation is useless, so skip it.
				num_factoring_failures++;
			}
			else if (other_poly != NULL && !Timed_complete_other_side(index, other_side_list))
			{
				num_other_side_failures++;
			}
			else
			{
				if (other_poly == NULL)
				{
					Get_saved_list(index, other_side_list);
				}
				const string& saved_list = other_side_list;

				//The factor list for the batch checked side: large primes first, then the other factors.
				checked_side_list.clear();

				for (int j = 0; j < num_large_primes; j++)
				{
					if (large_primes[j] >= (unsigned int)min_saved_prime_limit)
					{
						append_prime(checked_side_list, large_primes[j]);
					}
				}

				for (int j = 0; j < num_factors; j++)
				{
					if (relation_factor_list[j] >= (unsigned int)min_saved_prime_limit)
					{
						append_prime(checked_side_list, relation_factor_list[j]);
					}
				}

				//The relation is a,b:rational factors:algebraic factors
				sprintf(ab_text, "%lld,%lld:", a_vals[index], b_vals[index]);
				relation_text = ab_text;

				if (side == 0)
				{
					//rational factors were just found, the algebraic ones were saved.
					relation_text += checked_side_list;
					relation_text += ':';
					relation_text += saved_list;
				}
				else
				{
					//algebraic factors were just found, the rational ones were saved.
					relation_text += saved_list;
					relation_text += ':';
					relation_text += checked_side_list;
				}

				relation_text += '\n';

				//print the found relation
				fputs(relation_text.c_str(), outfile);

				num_good++;
			}
		}

		//reset for next iteration.
		*inputs[i] = 1;
		a_vals[index] = 0;
		b_vals[index] = 0;
		is_smooth[index] = false;
		lp1[index] = 0;
		lp2[index] = 0;
	}

	fclose(outfile);

	profile_add(PROF_WRITE, t_start);
	profile_cycles[PROF_WRITE] -= profile_cycles[PROF_FACTOR] + profile_cycles[PROF_OTHER_SIDE] - factor_start; //factoring is counted separately.

	if (three_lp_stage)
	{
		cout << "Ran a 3LP batch of " << entries << " candidates, found " << num_good << " relations" << endl;
	}
	else
	{
		cout << "Ran a batch, found " << num_good << " smooth numbers" << endl;
	}
	num_relations_found += num_good;
	entries = 0;
}

long long Batch_smooth::Get_num_relations_found()
{
	return num_relations_found + ((second_stage != NULL) ? second_stage->Get_num_relations_found() : 0);
}