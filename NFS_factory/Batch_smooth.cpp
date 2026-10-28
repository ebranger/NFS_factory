#include "Batch_smooth.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <cerrno>
#include <time.h>
#include "cofactorize.h"

#include "smintfact.h"
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

#define LP_PRIMALITY_REPS 10 //Miller-Rabin rounds used before accepting a cofactor as a large prime in a relation.

#define WRAPAROUND_MIN_LIMBS 512 //Use wraparound products in the remainder tree for nodes of at least this many limbs. Tuned on an i7-14700F.
#define FFT_MULMOD_2EXPP1_CUTOFF 256 //From the MPIR build (gmp-mparam.h). Below this, MPIR does not use the FFT for products mod B^n + 1.

static void append_prime(string& list, unsigned long long prime);

//MPIR functions for multiplication modulo B^n + 1 with the FFT. They are declared in gmp-impl.h, which is not part of the installed headers.
extern "C"
{
	int __gmpn_mulmod_Bexpp1_fft(mp_ptr op, mp_size_t pl, mp_srcptr n, mp_size_t nl, mp_srcptr m, mp_size_t ml);
	mpir_si __mpir_fft_adjust_limbs(mp_size_t limbs);
}

//Initialize the batch smoothness class
Batch_smooth::Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side, int print_limit, Polynomial* other_side_poly, int max_tree_levels)
{
	//init.
	temp_number = new mpz_class;

	long long max_prime_limit = 1LL << lbp;

	//if fblim and 2^lbp match, then we batch-check all numbers against the full primorial and do no cofactorization of any remaining numbers.
	// fblim > 2^lbp makes little sense? so assume equal, otherwise we will end up printing relations with some factors too large. Maybe check this in input file and warn user?
	if (max_prime_limit >= fblim)
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

	side = batch_checking_side;

	prime_product = new mpz_class();
	*prime_product = 1;

	Setup_prime_product(fblim);
	max_batch_prime = fblim;

	//For handling 2lp cofactorization. 
	one_lp_bit_limit = lbp;
	two_lp_min_limit = 2 * (long long)floor(log2((double)fblim));
	max_cofactor_bit_size = mfb;
	
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
		//Too big to be split into two primes both less than lp_limit
		return false;
	}

	if (cofactor_bits > two_lp_min_limit)
	{
		//Can potentially split into two large primes, small enough to be useful
		if (mpz_probab_prime_p(number.get_mpz_t(), 1))
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
			//either squfof failed or it is too big for squfof, so do qs.
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
		if (!mpz_probab_prime_p(factor1->get_mpz_t(), LP_PRIMALITY_REPS))
		{
			return false;
		}
		if (!mpz_probab_prime_p(factor2->get_mpz_t(), LP_PRIMALITY_REPS))
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
	if (number > max_batch_prime && mpz_probab_prime_p(number.get_mpz_t(), LP_PRIMALITY_REPS))
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
		num_factors = factor(relation_factor_list, number.get_mpz_t(), 0);

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

//Perform the remainder tree calculations.
//Core parts of this function was taken from the CADO-NFS package.
//Use a scaled remainder tree for speed.
void Batch_smooth::Do_batch_check()
{
	int limit;

	if (entries == 0)
	{
		//Nothing to check.
		return;
	}

	unsigned long long t = profile_now();

	//set up polynomial values
	for (int index = 0; index < entries; index++)
	{
		//use different poly evaluation if polynomial has form k*x^4+l which is faster.
		if (using_snfs_deg4_binomial)
		{
			poly->evaluate_binomial_polynomial_homogeneous_SNFS_deg4(inputs[start_entry + index], a_vals[index], b_vals[index], temp_number);
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

	Q = *prime_product % *inputs[0]; //do a mod since the prime product may well be much larger than the top number of the tree.
	t = profile_add(PROF_TOP_MOD, t);
	mpz_mul_2exp(Q.get_mpz_t(), Q.get_mpz_t(), product_bits[0] + guard);
	mpz_tdiv_q(inputs[0]->get_mpz_t(), Q.get_mpz_t(), inputs[0]->get_mpz_t()); // now the root of the tree is prime_product / product_of_input_numbers * 2^(bits in product of input number + a few guard bits)
	t = profile_add(PROF_TOP_SCALE, t);

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
				Middle_product(wrap_result[0].get_mpz_t(), inputs[(i - 1) / 2]->get_mpz_t(), inputs[i + 1]->get_mpz_t(), parent_bits - product_bits[i], parent_bits + guard);
				Middle_product(wrap_result[1].get_mpz_t(), inputs[(i - 1) / 2]->get_mpz_t(), inputs[i]->get_mpz_t(), parent_bits - product_bits[i + 1], parent_bits + guard);
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

	Output_solutions();

	delete[] product_bits;
	return;
}

//Compute the bits from shift up to precision of (y * x) mod 2^precision, i.e. (y * x mod 2^precision) >> shift.
//y must be less than 2^precision, and x less than 2^(shift + 1).
//The product is computed modulo B^n + 1, with n a little more than precision bits, using the FFT of MPIR. That is
//cheaper than the full product. The part of the product above B^n wraps around and is subtracted from the low bits,
//below the bits that are needed, so the result can be one too small. The guard bits of the scaled remainder tree
//absorb such errors, just as they absorb the truncation errors of the shifts.
void Batch_smooth::Middle_product(mpz_ptr result, mpz_srcptr y, mpz_srcptr x, mp_bitcnt_t shift, mp_bitcnt_t precision)
{
	mp_size_t limbs = (mp_size_t)((precision + 2 * GMP_LIMB_BITS - 1) / GMP_LIMB_BITS); //at least 64 bits more than precision.
	if (limbs >= FFT_MULMOD_2EXPP1_CUTOFF)
	{
		limbs = (mp_size_t)__mpir_fft_adjust_limbs(limbs);
	}

	if (wrap_buffer.size() < (size_t)limbs + 1)
	{
		wrap_buffer.resize(limbs + 1);
	}
	mp_ptr product = wrap_buffer.data();

	if (__gmpn_mulmod_Bexpp1_fft(product, limbs, y->_mp_d, abs(y->_mp_size), x->_mp_d, abs(x->_mp_size)) != 0)
	{
		//The product is B^n, i.e. -1 mod B^n + 1, which only happens if the needed bits are all zero.
		mpz_set_ui(result, 0);
		return;
	}

	//Shift down the needed bits.
	mp_bitcnt_t out_bits = precision - shift;
	mp_size_t out_limbs = (mp_size_t)((out_bits + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS);
	mp_size_t first_limb = (mp_size_t)(shift / GMP_LIMB_BITS);
	mp_size_t source_limbs = (mp_size_t)((precision + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS) - first_limb;

	if (result->_mp_alloc < source_limbs)
	{
		mpz_realloc2(result, source_limbs * GMP_LIMB_BITS);
	}

	if (shift % GMP_LIMB_BITS != 0)
	{
		mpn_rshift(result->_mp_d, product + first_limb, source_limbs, (unsigned int)(shift % GMP_LIMB_BITS));
	}
	else
	{
		mpn_copyi(result->_mp_d, product + first_limb, source_limbs);
	}

	if (out_bits % GMP_LIMB_BITS != 0)
	{
		result->_mp_d[out_limbs - 1] &= (((mp_limb_t)1) << (out_bits % GMP_LIMB_BITS)) - 1;
	}

	while (out_limbs > 0 && result->_mp_d[out_limbs - 1] == 0)
	{
		out_limbs--;
	}
	result->_mp_size = (int)out_limbs;
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
		int num_other_factors = factor_smooth(other_side_factor_list, norm.get_mpz_t(), MAX_FACTORS_IN_RELATION);
		if (num_other_factors < 0)
		{
			num_other_factors = factor(other_side_factor_list, norm.get_mpz_t(), 1);
		}

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

			//factor the relation. First, divide out already factored large primes.
			num_factors = 0;

			if (lp1[index] > 1)
			{
				mpz_divexact_ui(relation_value->get_mpz_t(), relation_value->get_mpz_t(), lp1[index]);
			}

			if (lp2[index] > 1)
			{
				mpz_divexact_ui(relation_value->get_mpz_t(), relation_value->get_mpz_t(), lp2[index]);
			}

			//run Pollard rho factoring on the remaining part of the relation. Since it is smooth Pollard rho should be a decent choice.
			if (*relation_value > 1)
			{
				unsigned long long tf = profile_now();
				num_factors = factor_smooth(relation_factor_list, relation_value->get_mpz_t(), MAX_FACTORS_IN_RELATION);
				if (num_factors < 0)
				{
					//Too large for the fast method, use the general one.
					num_factors = factor(relation_factor_list, relation_value->get_mpz_t(), 1);
				}
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

				if (lp1[index] > 1)
				{
					append_prime(checked_side_list, lp1[index]);
				}

				if (lp2[index] > 1)
				{
					append_prime(checked_side_list, lp2[index]);
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

	cout << "Ran a batch, found " << num_good << " smooth numbers" << endl;
	num_relations_found += num_good;
	entries = 0;
}

int Batch_smooth::Get_num_relations_found()
{
	return num_relations_found;
}