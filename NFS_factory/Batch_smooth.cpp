#include "Batch_smooth.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <cerrno>
#include <time.h>
#include "cofactorize.h"

#include "smintfact.h"

using namespace std;

#define MAX_TREE_LEVELS 20
//max 22 requires ~3GB of RAM. However, for SNFS-200 I seem to get the same speed for lower size, which saves a lot of RAM. 20 is ~7% slower, but saved 2GB of RAM, so worth it.
//For larger factorizations, consider increasing this a little, but for every +1 the memory requirement for the tree roughly doubles.

#define MAX_SQUARINGS_IN_CHECK 1 //Bernsteins algorithm does squaring at the bottom nodes, to handle powers of primes. This is the largest power checked.
//may want to do a small-prime check spearately, so that only the larger primes are still present when checking.
// Note: max squarings 2 means max power is 2^2=4.
//Currently: do only 1 squaring, it misses a few very rare relations but increase speed slightly.

#define MAX_CHARS_IN_RELATION_LIST 150 //from the relations i have, this should be enought space to hold any relation. 

//Initialize the batch smoothness class
Batch_smooth::Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side)
{
	//init.
	temp_number = new mpz_class;

	long long max_prime_limit = pow(2.0, lbp);

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

	//For handling 2lp cofactorization. 
	one_lp_bit_limit = lbp;
	two_lp_min_limit = 2 * floor(log2(fblim));
	max_cofactor_bit_size = mfb;
	
	entries = 0;

	//Estiamte a good tree size for the remainder tree
	tree_size = (int(log2(fblim))) - 5; // rough guess, works for deg-4 SNFS quite well.
	tree_size = min(MAX_TREE_LEVELS, tree_size); 

	start_entry = pow(2, (double)(tree_size - 1)) - 1;
	max_entries_before_check = pow(2, (double)tree_size - 1);
	save_filename = filename;
	num_relations_found = 0;

	prime_product = new mpz_class();
	*prime_product = 1;

	Setup_prime_product(fblim);
	product_set_up = true;
	max_batch_prime = fblim;

	Setup_memory(tree_size, max_bits_in_input);

	poly = input_poly;

	if (poly->poly_is_binomial_SNFS_deg4())
	{
		using_snfs_deg4_binomial = true;
	}
	else
	{
		using_snfs_deg4_binomial = false;
	}

	side = batch_checking_side;

	relation_factor_list = new unsigned int[20]; //I expect more like 10 factors maximum to be stored, but better not run out of space!

	precheck = true; //Debug, remove all tiny factors before batch-checking. Does not really affect the results much, though one or two more relations are found. And solves some crashes in the squfof code...

	factor1 = new mpz_class;
	factor2 = new mpz_class;
}

//Destructor, free all memory that was requested.
Batch_smooth::~Batch_smooth(void)
{
	Output_solutions();
	delete[] inputs;
	delete prime_product;
	delete no_print_product;
	delete predivide_tiny_product;

	for (int i = 0; i < num_relations_per_batch; i++)
	{
		delete list_of_factors[i];
	}
	delete[] list_of_factors;
	delete lp1;
	delete lp2;
	delete factor1;
	delete factor2;
}

//Setup the primorials to be used.
void Batch_smooth::Setup_prime_product(int smooth_bound)
{
	//So apparently mpir knows primorials. Problem solved!

	if (smooth_bound < 1000)
	{
		//I expect a minimum bound of ~2^20, anything less is probably a bad idea. 
		//Use 1000 as a basic sanity check. 
		return;
	}

	mpz_primorial_ui(prime_product->get_mpz_t(), smooth_bound);
	mpz_primorial_ui(temp_number->get_mpz_t(), int(sqrt((double)smooth_bound))); //Extra small primes. this saves having to do some squarings in the batch tree bottom node, and adds very few bits to the product.
	*prime_product = *prime_product * *temp_number;
	mpz_primorial_ui(temp_number->get_mpz_t(), int(sqrt(sqrt((double)smooth_bound)))); // Extra even smaller primes.
	*prime_product = *prime_product * (*temp_number * *temp_number);

	predivide_tiny_product = new mpz_class;
	*predivide_tiny_product = "10821610800"; // trial divide out the smallest primes that may occur to some higher power. This tests 2*2*2*2*3*3*3*5*5*7*7*11*11*13*13, for all 4-bit primes at least squared.

	no_print_product = new mpz_class;
	mpz_primorial_ui(no_print_product->get_mpz_t(), 1000); // For removing all factors less than 1000, since they will not be printed.

	long size = mpz_sizeinbase(prime_product->get_mpz_t(), 2);
	cout << "Size of prime product: " << size << " bits." << endl;

	//test for removing small factors before the larger multiplications happen.
	//small_factors = new mpz_class;
	//mpz_primorial_ui(small_factors->get_mpz_t(), 262144 * 4);

	return;
}

//allocate sufficient memory to hold the remainder trees. Since their size may change from run to run, allocate enough space that no realloc needs to be done. 
//to check: is the MPIR mpz_class memory mangement efficient enough that this can be skipped?
void Batch_smooth::Setup_memory(int tree_levels, int max_bits_in_input)
{
	//allocate enough memory to hold all the numbers in the product tree without having to increase any mpz size
	int temp = pow(2, (double)(tree_levels));
	inputs = new mpz_class*[temp - 1];
	lp1 = new unsigned long long[temp];
	lp2 = new unsigned long long[temp];
	int limit;

	temp = temp / 2;
	num_relations_per_batch = temp;
	a_vals = new long long[num_relations_per_batch];
	b_vals = new long long[num_relations_per_batch];

	list_of_factors = new char*[num_relations_per_batch];

	for (int i = 0; i < temp; i++)
	{
		//allocate space for storing the relation factor list
		list_of_factors[i] = new char[MAX_CHARS_IN_RELATION_LIST];

		//set large prime list to 0
		lp1[i] = 0;
		lp2[i] = 0;
	}

	temp = max_bits_in_input;

	for (int level = tree_levels; level >= 1; level--)
	{
		limit = pow(2, (double)(level - 1)) - 1;

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

//Add one relation to be batch-smoothness checked.
void Batch_smooth::Add_number(long long a_val, long long b_val, char* relation_factor_list)
{
	//add one relation for batch-factoring.

	a_vals[entries] = a_val;
	b_vals[entries] = b_val;
	strcpy(list_of_factors[entries], relation_factor_list);

	entries++;

	//Is the tree full?
	if (entries == max_entries_before_check) //entries points at first empty position.
	{
		Do_batch_check();
	}
}

//Check whether a remainder factors into two numbers less than the lp limit. 
//If it does, it is a valid relation!
bool Batch_smooth::check_two_lp_smooth(mpz_class number, unsigned int index)
{
	int cofactor_bits = mpz_sizeinbase(number.get_mpz_t(), 2);

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

		//Update: use SIQS, so 0 means fail and 1 means sucess, the name num_factors is now misleading...
		num_factors = tinyqs(number.get_mpz_t(), factor1->get_mpz_t(), factor2->get_mpz_t());

		if (num_factors == 0)
		{
			//unsuccesful
			return false;
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
		if (mpz_sizeinbase(factor1->get_mpz_t(), 2) > one_lp_bit_limit)
		{
			return false;
		}
		if (mpz_sizeinbase(factor2->get_mpz_t(), 2) > one_lp_bit_limit)
		{
			return false;
		}

		//save the two found factors. I do not expect three factors to ever happen, though that may cause problems here...
		lp1[index] = factor1->get_ui();
		lp2[index] = factor2->get_ui();

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
	if (mpz_probab_prime_p(number.get_mpz_t(), 1) && number > max_batch_prime)
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

		for (int i = 0; i <= num_factors; i++)
		{
			if (log2(relation_factor_list[i]) >= one_lp_bit_limit)
			{
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
	int index;

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

		if (precheck == true)
		{
			//divide out tiny factors that can easily occur to high powers.

			//mpz_gcd(temp_number->get_mpz_t(), inputs[start_entry + index]->get_mpz_t(), predivide_tiny_product->get_mpz_t());
			mpz_gcd(temp_number->get_mpz_t(), inputs[start_entry + index]->get_mpz_t(), predivide_tiny_product->get_mpz_t());
			//cout << "begin: \n" <<  *inputs[start_entry + entries] << "\n" << *small_prime_product << "\n" << temp_number << endl;
			while (*temp_number > 1)
			{
				mpz_divexact(inputs[start_entry + index]->get_mpz_t(), inputs[start_entry + index]->get_mpz_t(), temp_number->get_mpz_t());
				//*inputs[start_entry + entries] = *inputs[start_entry + entries] / *temp_number;
				mpz_gcd(temp_number->get_mpz_t(), inputs[start_entry + index]->get_mpz_t(), temp_number->get_mpz_t());
				//cout << *inputs[start_entry + entries] << endl;
			}
		}
	}

	//set up product tree.
	for (int level = tree_size; level > 1; level--)
	{
		limit = pow(2, (double)(level - 2)) - 1;

		for (int i = limit; i <= 2 * limit; i++)
		{
			*inputs[i] = *inputs[2 * i + 1] * *inputs[2 * i + 2];
		}
	}

	long size = mpz_sizeinbase(inputs[0]->get_mpz_t(), 2);
	cout << "Size of top node of product tree: " << size << " bits." << endl;


	//The core remainder tree calculations are based on the CADO-nfs batch function

	long h = tree_size + 2; // extra fudge bits.
	unsigned long i, j, guard, *product_bits;
	product_bits = new unsigned long[2 * num_relations_per_batch];

	guard = h;

	//Get the size of all products, they are needed for shifts later.
	for (i = 0; i < 2 * num_relations_per_batch - 1; i++)
	{
		product_bits[i] = mpz_sizeinbase(inputs[i]->get_mpz_t(), 2);
	}

	mpz_class Q;

	Q = *prime_product % *inputs[0]; //do a mod since the prime product may well be much larger than the top number of the tree.
	mpz_mul_2exp(Q.get_mpz_t(), Q.get_mpz_t(), product_bits[0] + guard);
	mpz_tdiv_q(inputs[0]->get_mpz_t(), Q.get_mpz_t(), inputs[0]->get_mpz_t()); // now the root of the tree is prime_product / product_of_input_numbers * 2^(bits in product of input number + a few guard bits)

	for (int level = 2; level <= tree_size - 1; level++)
	{
		limit = pow(2, (double)(level - 1)) - 1;

		for (int i = limit; i <= 2 * limit; i = i + 2)
		{
			//calculate the next level of the tree. 
			// so if the parent node is Num / (leaf1 * leaf2), multiply by the leaves and right-shift since the fraction is now smaller, requiring fewer bits of storage.

			//Q = *inputs[i];
			//*inputs[i] = *inputs[(i - 1) / 2] * *inputs[i + 1];
			//*inputs[i + 1] = *inputs[(i - 1) / 2] * Q;

			*inputs[i] = *inputs[(i - 1) / 2] * *inputs[i];
			*inputs[i + 1] = *inputs[(i - 1) / 2] * *inputs[i + 1];
			mpz_swap(inputs[i]->get_mpz_t(), inputs[i + 1]->get_mpz_t());

			mpz_tdiv_r_2exp(inputs[i]->get_mpz_t(), inputs[i]->get_mpz_t(), product_bits[(i - 1) / 2] + guard);
			mpz_div_2exp(inputs[i]->get_mpz_t(), inputs[i]->get_mpz_t(), product_bits[(i - 1) / 2] - product_bits[i]);

			mpz_tdiv_r_2exp(inputs[i + 1]->get_mpz_t(), inputs[i + 1]->get_mpz_t(), product_bits[(i - 1) / 2] + guard);
			mpz_div_2exp(inputs[i + 1]->get_mpz_t(), inputs[i + 1]->get_mpz_t(), product_bits[(i - 1) / 2] - product_bits[i + 1]);
		}
	}

	//for the bottom node, i.e. the relation numbers, store the Prime_product % number in Q, to keep the number (in input_number[]). Then it can be used later for squaring tests mod the number.
	for (i = num_relations_per_batch - 1; i < 2 * num_relations_per_batch - 1; i = i + 2)
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
			if (!check_two_lp_smooth(Q, i - num_relations_per_batch + 1))
			{
				//not smooth.
				a_vals[i - num_relations_per_batch + 1] = 0;
				b_vals[i - num_relations_per_batch + 1] = 0;
			}
		}
		else
		{
			if (Q > 0)
			{
				//Do not do 2lp testing. If Q = 0 then the number is smooth

				//not smooth.
				a_vals[i - num_relations_per_batch + 1] = 0;
				b_vals[i - num_relations_per_batch + 1] = 0;
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
			if (!check_two_lp_smooth(Q, i - num_relations_per_batch + 2))
			{
				//not smooth.
				a_vals[i - num_relations_per_batch + 2] = 0;
				b_vals[i - num_relations_per_batch + 2] = 0;
			}
		}
		else
		{
			if (Q > 0)
			{
				//Do not do 2lp testing. If Q = 0 then the relation is smooth

				//not smooth.
				a_vals[i - num_relations_per_batch + 2] = 0;
				b_vals[i - num_relations_per_batch + 2] = 0;
			}
		}
	}

	Output_solutions();

	delete product_bits;
	return;
}

//For each relation passing the batch smoothnes checking, complete the relation and print it.
void Batch_smooth::Output_solutions()
{
	//Print all found relations to file
	int limit;
	limit = pow(2, (double)(tree_size - 1)) - 1;
	int index;

	char* relation_text = new char[200]; // i expect that a relation fits in a maximum of 130 bytes, but would rather not run out of space.
	char* text_ptr;
	int char_written;

	int num_good = 0;

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
		if (a_vals[index] != 0 && index < entries)
		{
			relation_value = inputs[i];

			//negative values are handled by msieve, ignore here. 
			if (*relation_value < 0)
			{
				*relation_value = *relation_value * -1;
			}

			//divide out factors < 1000 that will not be printed. Msieve reconstructs these.
			mpz_gcd(temp_number->get_mpz_t(), relation_value->get_mpz_t(), no_print_product->get_mpz_t());

			while (*temp_number > 1)
			{
				mpz_divexact(relation_value->get_mpz_t(), relation_value->get_mpz_t(), temp_number->get_mpz_t());
				//*relation_value = *relation_value / *temp_number;
				mpz_gcd(temp_number->get_mpz_t(), relation_value->get_mpz_t(), temp_number->get_mpz_t());
			}

			//Now to factor the relation value
			if (*relation_value > 1)
			{
				//factor the relation. First, check for already factored 2lp primes.
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
				num_factors += factor(relation_factor_list + num_factors, relation_value->get_mpz_t(), 1);
			}
			else
			{
				//very smooth relation, all factors less than 1000. No factors to print.
				num_factors = 0;
			}

			if (num_factors == -1)
			{
				cout << "failed to factor a relation!" << endl;
				cout << a_vals[index] << " " << b_vals[index] << endl;
			}

			//Print the relation to file.
			text_ptr = relation_text;

			char_written = sprintf(text_ptr, "%lld", a_vals[index]);

			text_ptr += char_written;
			*text_ptr = ',';
			text_ptr++;

			char_written = sprintf(text_ptr, "%lld", b_vals[index]);

			text_ptr += char_written;
			*text_ptr = ':';
			text_ptr++;
			*text_ptr = '\0';

			if (side == 0)
			{
				//rational relations to be printed first. These should just have been found and stored in the list. 
			}
			else
			{
				//algebraic relations are found, so print the saved rational factors first.
				strcpy(text_ptr, list_of_factors[index]);
				text_ptr += strlen(text_ptr);
				*text_ptr = ':';
				text_ptr++;
				*text_ptr = '\0';
			}

			//Handle large primes that may have been saved earlier.
			if (lp1[index] > 1)
			{
				char_written = sprintf(text_ptr, "%X", lp1[index]);
				text_ptr += char_written;
				*text_ptr = ',';
				text_ptr++;
				lp1[index] = 0;
			}

			if (lp2[index] > 1)
			{
				char_written = sprintf(text_ptr, "%X", lp2[index]);
				text_ptr += char_written;
				*text_ptr = ',';
				text_ptr++;
				lp2[index] = 0;
			}

			//Print the found factors for the smoothnes.checked and factored relation value.
			for (int i = 0; i < num_factors; i++)
			{
				char_written = sprintf(text_ptr, "%X", relation_factor_list[i]);

				text_ptr += char_written;
				if (i < num_factors)
				{
					*text_ptr = ',';
					text_ptr++;
				}
			}

			//remove the final comma that was printed.
			text_ptr--;
			*text_ptr = '\0';

			if (side == 0)
			{
				*text_ptr = ':';
				text_ptr++;

				//rational relations were just found and printed, now fill in the saved algebraic ones.
				//TODO: test this program when batch checking rational values, and test the outputting done here
				strcpy(text_ptr, list_of_factors[index]);
				text_ptr += strlen(text_ptr);
			}
			else
			{
				//algebraic relations are found, so now all rational and algebraic relations have already been printed.

			}

			*text_ptr = '\n';
			text_ptr++;
			*text_ptr = '\0';

			//print the found relation
			fputs(relation_text, outfile);

			num_good++;
		}

		//reset for next iteration.
		*inputs[i] = 1;
		a_vals[index] = 0;
		b_vals[index] = 0;

	}

	fclose(outfile);

	cout << "Ran a batch, found " << num_good << " smooth numbers" << endl;
	num_relations_found += num_good;
	entries = 0;
	delete relation_text;
}

int Batch_smooth::Get_num_relations_found()
{
	return num_relations_found;
}