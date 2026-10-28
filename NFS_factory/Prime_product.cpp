#include "Prime_product.h"
#include <vector>
#include <iostream>
using namespace std;

//Adapted from code used to generate a factorbase for a lattice siever, 

#define ERATOSTHENES_LINEWIDTH 16384 //Fits in cache, works for now.

void multiply_recursively(vector<uint32> primes, uint32 lowest_prime_index, uint32 highest_prime_index, mpz_class* product, uint32 square_lim, uint32 double_square_lim)
{
	if (highest_prime_index - lowest_prime_index >= 10)
	{
		int mid_index = (highest_prime_index + lowest_prime_index) / 2;
		multiply_recursively(primes, lowest_prime_index, mid_index, product, square_lim, double_square_lim);
		multiply_recursively(primes, mid_index + 1, highest_prime_index, product, square_lim, double_square_lim);
	}
	else 
	{
		for (int i = lowest_prime_index; i <= highest_prime_index; i++)
		{
			*product = *product * primes[i];
			if (primes[i] >= square_lim)
			{
				*product = *product * primes[i];
			}
			if (primes[i] >= double_square_lim)
			{
				*product = *product * primes[i];
				*product = *product * primes[i];
			}
		}
	}
}

void setup_prime_product(Polynomial* poly, uint32 maxprime, mpz_class* prime_product)
{
	*prime_product = 1;
	int num_primes = 0;

	if (maxprime < 10)
	{
		cout << "Error: Prime_product was called with a too small maximum prime!" << endl;
		return;
	}

	if (poly->get_degree() <= 1)
	{
		mpz_class* temp_number = new mpz_class;
		*temp_number = 1;
		//Degree 1 polynomial, all primes can divide the norm. Just do a primorial.
		mpz_primorial_ui(prime_product->get_mpz_t(), maxprime);
		mpz_primorial_ui(temp_number->get_mpz_t(), int(sqrt((double)maxprime))); //Extra small primes. this saves having to do some squarings in the batch tree bottom node, and adds very few bits to the product.
		*prime_product = *prime_product * *temp_number;
		mpz_primorial_ui(temp_number->get_mpz_t(), int(sqrt(sqrt((double)maxprime)))); // Extra even smaller primes.
		*prime_product = *prime_product * (*temp_number * *temp_number);
		delete temp_number;
		return;
	}

	//Degree 2 or higher polynomial, need to find which primes can divide the norm.

	uint32 max_sieve_prime = int(sqrt((double)maxprime)) + 1;
	uint32 number_of_sieve_lines = uint32(ceil((double)maxprime / (ERATOSTHENES_LINEWIDTH * 2)));

	vector<uint32> sieve_primes;
	vector<uint32> primelist;

	//Skip all entries divisible by 2. the line is then the numbers 1,3,5,7....
	uint8* sieve_line = new uint8[ERATOSTHENES_LINEWIDTH];

	for (uint32 i = 0; i < ERATOSTHENES_LINEWIDTH; i++)
	{
		sieve_line[i] = 1;
	}

	sieve_line[0] = 0; // 1 is no prime
	uint32 current_prime;
	uint32 sieve_location;

	uint32* number_of_roots = new uint32;
	uint32* rootlist = new uint32[poly->get_degree()]; //no more than this number of roots can occur

	//poly->get_polyroots_modp(number_of_roots, rootlist, 2); //handle 2 separately, since later only odd numbers are sieved
	*number_of_roots = 0;

	if (*number_of_roots > 0)
	{
		primelist.push_back(2);
		num_primes++;
	}

	//handle the first line separately. expect that all primes < sqrt(maxprime) occur in this line, these primes will be sieved later.
	for (uint32 line_index = 1; line_index < ERATOSTHENES_LINEWIDTH; line_index++)
	{
		if (sieve_line[line_index] == 1)
		{
			current_prime = line_index * 2 + 1;

			if (current_prime < max_sieve_prime)
			{
				//save small prime for later, it will be used for sieving.
				sieve_primes.push_back(current_prime);
			}
			sieve_location = line_index + current_prime;

			//do sieving! 
			while (sieve_location < ERATOSTHENES_LINEWIDTH)
			{
				sieve_line[sieve_location] = 0;
				sieve_location += current_prime;
			}

			//append the prime and roots to the factor base.
			poly->get_polyroots_modp(number_of_roots, rootlist, current_prime);

			if (*number_of_roots > 0)
			{
				primelist.push_back(current_prime);
				num_primes++;
			}
		}
	}

	uint32 sieving_prime = 0;

	//sieve the remaining lines 
	for (uint32 line = 1; line < number_of_sieve_lines; line++)
	{
		//reset
		for (uint32 i = 0; i < ERATOSTHENES_LINEWIDTH; i++)
		{
			sieve_line[i] = 1;
		}

		for (uint32 prime_index = 0; prime_index < sieve_primes.size(); prime_index++)
		{
			sieving_prime = sieve_primes[prime_index];
			//The first location that is divisible by the prime
			sieve_location = ((sieving_prime - 1) / 2 + (sieving_prime - (ERATOSTHENES_LINEWIDTH % sieving_prime)) * line) % sieving_prime;

			while (sieve_location < ERATOSTHENES_LINEWIDTH)
			{
				sieve_line[sieve_location] = 0;
				sieve_location += sieving_prime;
			}
		}

		//find the primes
		for (uint32 i = 0; i < ERATOSTHENES_LINEWIDTH; i++)
		{
			current_prime = 2 * (i + ERATOSTHENES_LINEWIDTH * line) + 1;

			if (sieve_line[i] == 1 && current_prime < maxprime)
			{
				//Found a prime that may be in the fb.
				poly->get_polyroots_modp(number_of_roots, rootlist, current_prime);

				if (*number_of_roots > 0)
				{
					primelist.push_back(current_prime);
					num_primes++;
				}
			}
		}
	}

	multiply_recursively(primelist, 0, num_primes, prime_product, int(sqrt(maxprime)), int(sqrt(sqrt(maxprime))));
	
	//Cleanup
	delete number_of_roots;
	delete rootlist;
	delete sieve_line;
	primelist.empty();
}