/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Erik Branger. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may
benefit from your work.

Some parts of the code (and also this header), included in this
distribution have been reused from other sources. In particular I
have benefitted greatly from the work of Jason Papadopoulos's msieve,
the batch factorization code in CADO-NFS, and relation factoring code
found in gnfs-lasieve by Franke and Kleinjung.
----------------------------------------------------------------------*/


#include <iostream>
#include <fstream>
#include "mpirxx.h"
#include <math.h>
#include <string.h>
#include "Polynomial.h"
#include "Batch_smooth.h"
#include <time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>

using namespace std;

#define MAX_LINE_LENGTH 1024 //Longer lines in a relation file are skipped. gnfs-lasieve relations are typically less than 200 characters.
#define MAX_BINARY_PRIMES 32 //Maximum number of primes in one relation in a binary relation file. More than this is treated as a corrupt file.

struct nfs_parameters {
	//only those possibly needed. In principle only one side is needed...
	Polynomial* rational_polynomial, *algebraic_polynomial;
	int lbpr, lbpa;
	int smoothness_checking_side; // 0 for rational, 1 for algebraic.
	long long maxa, maxb; // Relations with larger |a| or b are skipped.
	long long maxrels; // Stop after reading this many relations. 0 means no limit.
	long long rlim, alim;
	int mfbr, mfba;
	int max_rational_coefficient_bits;
	int max_algebraic_coefficient_bits;
	int print_limit; // Skip printing all primes less than this.
	bool cado_format; // Print all primes on both sides, with multiplicity, as CADO-NFS requires.
};

//Remove leading and trailing whitespace.
string trim(const string& text)
{
	size_t first = text.find_first_not_of(" \t\r\n");
	if (first == string::npos)
	{
		return "";
	}
	size_t last = text.find_last_not_of(" \t\r\n");
	return text.substr(first, last - first + 1);
}

//Function for reading a polynomial. Read only the parameters of interest to this program and store in an nfs_parameter struct.
//Returns false if the file could not be read, or if a value could not be parsed.
bool read_polynomial(string filename, nfs_parameters* params)
{
	ifstream read;
	read.open(filename, ifstream::in);

	if (read.fail())
	{
		cout << "Failed to open input polynomial file." << endl;
		return false;
	}

	//values indicating nothing read.
	params->lbpa = 0;
	params->lbpr = 0;
	params->mfba = 0;
	params->mfbr = 0;
	params->alim = 0;
	params->rlim = 0;
	params->smoothness_checking_side = 2;

	params->max_algebraic_coefficient_bits = 0;
	params->max_rational_coefficient_bits = 0;

	string line;
	string data;
	mpz_class number;
	int line_number = 0;

	while (getline(read, line))
	{
		line_number++;
		line = trim(line);

		if (line.length() == 0 || line[0] == '#')
		{
			continue;
		}

		size_t position = line.find(':');
		if (position == string::npos)
		{
			continue;
		}
		data = trim(line.substr(0, position));
		line = trim(line.substr(position + 1));

		try
		{
			//For rational root given as x - m in the input. This is a short form occasionally used.
			//Otherwise, both Y0 and Y1 needs to be set.
			if (data == "m")
			{
				number = line;
				number = number * -1;
				params->rational_polynomial->set_coeff(0, &number);
				params->rational_polynomial->set_coeff(1, 1);
				params->max_rational_coefficient_bits = (int)mpz_sizeinbase(number.get_mpz_t(), 2);
			}
			//Polynomials in GGNFS format.
			else if (data.length() == 2 && data[0] == 'c' && isdigit(data[1]))
			{
				number = line;
				params->algebraic_polynomial->set_coeff(data[1] - '0', &number);
				if ((int)mpz_sizeinbase(number.get_mpz_t(), 2) > params->max_algebraic_coefficient_bits)
				{
					params->max_algebraic_coefficient_bits = (int)mpz_sizeinbase(number.get_mpz_t(), 2);
				}
			}
			else if (data.length() == 2 && data[0] == 'Y' && isdigit(data[1]))
			{
				number = line;
				params->rational_polynomial->set_coeff(data[1] - '0', &number);
				if ((int)mpz_sizeinbase(number.get_mpz_t(), 2) > params->max_rational_coefficient_bits)
				{
					params->max_rational_coefficient_bits = (int)mpz_sizeinbase(number.get_mpz_t(), 2);
				}
			}
			//Factorbase limits on the algebraic side
			else if (data == "alim")
			{
				params->alim = stoll(line);
			}
			//Factorbase limits on the rational side
			else if (data == "rlim")
			{
				params->rlim = stoll(line);
			}
			//Maximum bits for large primes on the algebraic side
			else if (data == "lpba")
			{
				params->lbpa = stoi(line);
			}
			//Maximum bits for large primes on the rational side
			else if (data == "lpbr")
			{
				params->lbpr = stoi(line);
			}
			//Maximum bits for the remainder on the algebraic side after batch smoothness checking, to be factored by QS/squfof
			else if (data == "mfba")
			{
				params->mfba = stoi(line);
			}
			//Maximum bits for the remainder on the rational side after batch smoothness checking, to be factored by QS/squfof
			else if (data == "mfbr")
			{
				params->mfbr = stoi(line);
			}
			//The side to perform batch smoothnes checking on. 0 is rational, 1 is algebraic
			else if (data == "side")
			{
				params->smoothness_checking_side = stoi(line);
			}
			//maximum a value in a relation. Relations with larger (absolute) values are discarded
			else if (data == "maxa")
			{
				params->maxa = stoll(line);
			}
			//Maximum b value in a relation. Relations with larger values are discarded
			else if (data == "maxb")
			{
				params->maxb = stoll(line);
			}
			//Stop after reading this many relations from the relation file. Useful for testing.
			else if (data == "maxrels")
			{
				params->maxrels = stoll(line);
			}
			//skip printing small primes.
			else if (data == "print")
			{
				params->print_limit = stoi(line);
			}
			//Output format: msieve (default) or cado.
			else if (data == "format")
			{
				if (line == "cado")
				{
					params->cado_format = true;
				}
				else if (line == "msieve")
				{
					params->cado_format = false;
				}
				else
				{
					cout << "Unknown output format " << line << " on line " << line_number << " of the polynomial file. Use msieve or cado." << endl;
					return false;
				}
			}
		}
		catch (const exception&)
		{
			cout << "Could not read the value of " << data << " on line " << line_number << " of the polynomial file." << endl;
			return false;
		}
	}

	read.close();
	return true;
}

//Parse a decimal integer with an optional leading minus sign. Returns a pointer to the first character after the number,
//or NULL if there was no number, or it was too long to be a reasonable a or b value.
const char* parse_integer(const char* text, long long* value)
{
	bool negative = false;
	int digits = 0;

	if (*text == '-')
	{
		negative = true;
		text++;
	}

	*value = 0;
	while (isdigit(*text))
	{
		*value = 10 * *value + (*text - '0');
		text++;
		digits++;
	}

	if (digits == 0 || digits > 18)
	{
		return NULL;
	}

	if (negative)
	{
		*value = -*value;
	}
	return text;
}

//Check that a factor list only contains hexadecimal numbers separated by commas. An empty list is allowed.
bool valid_factor_list(const char* list)
{
	bool expect_digit = true;

	for (const char* c = list; *c != '\0'; c++)
	{
		if (isxdigit(*c))
		{
			expect_digit = false;
		}
		else if (*c == ',' && !expect_digit)
		{
			expect_digit = true;
		}
		else
		{
			return false;
		}
	}

	//No trailing comma, unless the list is empty.
	return !expect_digit || *list == '\0';
}

//Check if a relation is within the allowed range of a and b values.
bool relation_in_range(long long a, long long b, nfs_parameters* param)
{
	return llabs(a) <= param->maxa && b <= param->maxb;
}

//Function for reading relations to be processed.
int read_batch_from_file(string filename, Batch_smooth* batch, nfs_parameters* param)
{
	FILE* infile;

	infile = fopen(filename.c_str(), "r");

	if (infile == NULL)
	{
		cout << "Could not open input file." << endl;
		return -1;
	}

	char buffer[MAX_LINE_LENGTH];
	long long a, b;
	const char *tmp;
	char* ratlist, *alglist;

	long long relation_num = 0;
	long long failed_relations = 0;
	long long skipped_relations = 0;
	bool incomplete_last_line = false;

	//fgets and char checking is so much faster than stringstream...
	while (fgets(buffer, sizeof(buffer), infile) != NULL)
	{
		//read a,b,rational factor list and algebraic factor list from each line of the input.
		size_t length = strlen(buffer);

		if (length == 0 || buffer[length - 1] != '\n')
		{
			if (feof(infile))
			{
				//The last line has no newline. A siever that is stopped may leave such an incomplete line, so skip it.
				incomplete_last_line = true;
				break;
			}

			//Line too long. Skip the rest of it.
			int c;
			while ((c = fgetc(infile)) != EOF && c != '\n');
			failed_relations++;
			continue;
		}

		//Remove the line ending.
		buffer[--length] = '\0';
		if (length > 0 && buffer[length - 1] == '\r')
		{
			buffer[--length] = '\0';
		}

		if (length == 0 || buffer[0] == '#')
		{
			//Empty line or comment, so skip this line.
			continue;
		}

		//The line should be a,b:rational factor list:algebraic factor list
		ratlist = strchr(buffer, ':');
		alglist = (ratlist != NULL) ? strchr(ratlist + 1, ':') : NULL;

		if (alglist == NULL || strchr(alglist + 1, ':') != NULL)
		{
			failed_relations++;
			continue;
		}

		*ratlist++ = '\0';
		*alglist++ = '\0';

		//a is stored in decimal and can be negative, b is stored in decimal and is positive.
		tmp = parse_integer(buffer, &a);
		if (tmp == NULL || *tmp != ',')
		{
			failed_relations++;
			continue;
		}

		tmp = parse_integer(tmp + 1, &b);
		if (tmp == NULL || *tmp != '\0' || b <= 0)
		{
			failed_relations++;
			continue;
		}

		if (!valid_factor_list(ratlist) || !valid_factor_list(alglist))
		{
			failed_relations++;
			continue;
		}

		if (!relation_in_range(a, b, param))
		{
			skipped_relations++;
			continue;
		}

		//Add the relation to be batch smoothness checked. The factor list of the other side is saved for printing.
		if (!batch->Add_number(a, b, (batch->side == 0) ? alglist : ratlist))
		{
			failed_relations++;
			continue;
		}

		relation_num++;

		if (param->maxrels > 0 && relation_num >= param->maxrels)
		{
			cout << "Stopping after " << relation_num << " relations (maxrels)." << endl;
			break;
		}
	}

	fclose(infile);

	//Finish any remaining relations, even if the tree is not filled.
	batch->Do_batch_check();

	if (incomplete_last_line)
	{
		cout << "Skipped an incomplete last line in the relation file." << endl;
	}

	if (failed_relations > 0)
	{
		cout << "Failed to parse " << failed_relations << " relations out of a total of " << relation_num + failed_relations + skipped_relations << " relations." << endl;
	}

	if (skipped_relations > 0)
	{
		cout << "Skipped " << skipped_relations << " relations with a or b larger than maxa or maxb." << endl;
	}

	return batch->Get_num_relations_found();
}

//Convert a list of primes to a comma separated list of hexadecimal values for printing.
void get_relation_factor_list(unsigned int* primes, unsigned int num_primes, char* relation_line)
{
	char* text_ptr = relation_line;

	for (unsigned int i = 0; i < num_primes; i++)
	{
		if (i > 0)
		{
			*text_ptr++ = ',';
		}
		text_ptr += sprintf(text_ptr, "%X", primes[i]);
	}

	*text_ptr = '\0';
}

//Read relations from a binary format file.
int read_batch_from_binary_file(string filename, Batch_smooth* batch, nfs_parameters* param)
{
	//The binary format is (int)a, (uint)b, (uint)num_primes, (uint)prime1 ... (uint)prime_num_primes
	int a;
	unsigned int b;
	unsigned int num_primes;
	unsigned int primes[MAX_BINARY_PRIMES];
	char relation_line[MAX_BINARY_PRIMES * 9 + 1]; //At most 8 hex digits and a comma per prime.

	long long relation_num = 0;
	long long skipped_relations = 0;
	bool corrupt = false;

	ifstream infile;
	infile.open(filename, ios::binary);

	if (infile.fail())
	{
		cout << "Failed to open input saved relations file." << endl;
		return -1;
	}

	while (infile.peek() != EOF)
	{
		//read (a,b) value of relation, and how many primes are stored in the relation.
		infile.read(reinterpret_cast<char*>(&a), sizeof(int));
		infile.read(reinterpret_cast<char*>(&b), sizeof(unsigned int));
		infile.read(reinterpret_cast<char*>(&num_primes), sizeof(unsigned int));

		if (!infile || num_primes > MAX_BINARY_PRIMES)
		{
			corrupt = true;
			break;
		}

		//read the primes.
		infile.read(reinterpret_cast<char*>(primes), num_primes * sizeof(unsigned int));

		if (!infile)
		{
			corrupt = true;
			break;
		}

		if (!relation_in_range(a, b, param))
		{
			skipped_relations++;
			continue;
		}

		get_relation_factor_list(primes, num_primes, relation_line);

		//Add the relation to be batch smoothness checked.
		batch->Add_number(a, b, relation_line);
		relation_num++;

		if (param->maxrels > 0 && relation_num >= param->maxrels)
		{
			cout << "Stopping after " << relation_num << " relations (maxrels)." << endl;
			break;
		}
	}
	infile.close();

	//Check final relations, even if the tree is not filled.
	batch->Do_batch_check();

	if (corrupt)
	{
		cout << "The binary relation file is truncated or corrupt after " << relation_num + skipped_relations << " relations. The rest of the file was skipped." << endl;
	}

	if (skipped_relations > 0)
	{
		cout << "Skipped " << skipped_relations << " relations with a or b larger than maxa or maxb." << endl;
	}

	return batch->Get_num_relations_found();
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		cout << "Usage: NFS_factory <polynomial file name> <saved relation file name> <output file name>" << endl;
		return 0;
	}

	string inputfile, outputfile, polyfile;
	polyfile = argv[1];
	inputfile = argv[2];
	outputfile = argv[3];

	Batch_smooth* batch;

	clock_t time;

	nfs_parameters* param = new nfs_parameters();
	param->maxa = LLONG_MAX; //No limit unless set in the polynomial file.
	param->maxb = LLONG_MAX;
	param->maxrels = 0;

	param->print_limit = 1000; //Msieve default.
	param->cado_format = false;

	param->algebraic_polynomial = new Polynomial();
	param->rational_polynomial = new Polynomial();

	if (!read_polynomial(polyfile, param))
	{
		return -1;
	}

	//Check the polynomial to ensure input is reasonable.
	if (param->smoothness_checking_side != 0 && param->smoothness_checking_side != 1)
	{
		cout << "Could not find which side to check for smoothness. Set side to 0 (rational) or 1 (algebraic) in the polynomial file." << endl;
		return -1;
	}

	if (param->algebraic_polynomial->get_degree() == 0 || param->rational_polynomial->get_degree() == 0)
	{
		cout << "Failed to read the polynomials from the polynomial file." << endl;
		return -1;
	}

	//Values for the side to batch check.
	bool rational = (param->smoothness_checking_side == 0);
	long long lim = rational ? param->rlim : param->alim;
	int lbp = rational ? param->lbpr : param->lbpa;
	int mfb = rational ? param->mfbr : param->mfba;
	string side_name = rational ? "r" : "a";

	if (lim < 1000)
	{
		cout << "The factor base limit " << side_name << "lim is missing or too small in the polynomial file." << endl;
		return -1;
	}

	//Primes of 32 bits is the maximum that can be handled by the batch factoring code. Several changes are needed before 33 or larger will work (mainly in how an identified relation is fully factored).
	if (lim > 4294967295)
	{
		cout << "Batch smoothness checking can only be done for primes of 32 bits or less. Use a lower " << side_name << "lim" << endl;
		return -1;
	}

	if (lbp < 1 || lbp > 32)
	{
		cout << "The large prime bound lpb" << side_name << " is missing, or larger than 32 bits, in the polynomial file." << endl;
		return -1;
	}

	if (mfb < 1)
	{
		cout << "mfb" << side_name << " is missing from the polynomial file." << endl;
		return -1;
	}

	//Can only handle cofactorization with 2lp, so no point in checking numbers that are too big for this.
	if (lbp * 2 < mfb)
	{
		cout << "Three large prime cofactorization not supported. Use a lower mfb" << side_name << "." << endl;
		return -1;
	}

	//Estimate size needed to hold the value of the algebraic/rational polynomial values.
	int tree_bits;
	if (rational)
	{
		tree_bits = param->max_rational_coefficient_bits+30; //Assume max 30 bit a, b values in relation, when calculating the optimal tree size.
	}
	else
	{
		tree_bits = param->max_algebraic_coefficient_bits  + 30 * param->algebraic_polynomial->get_degree(); // Assume maxumum a, b is 30 bits and maximum coefficient is 30 bits (at bigger coefficients, batch factoring becomes too slow, and regular SNFS may be preferred? )
	}

	Polynomial* checked_poly = rational ? param->rational_polynomial : param->algebraic_polynomial;
	Polynomial* other_poly = rational ? param->algebraic_polynomial : param->rational_polynomial;

	if (param->cado_format)
	{
		//CADO-NFS needs all primes on both sides, so print everything and rebuild the factor list of the other side.
		param->print_limit = 1;
		cout << "Writing relations in CADO-NFS format, with all primes on both sides." << endl;
	}
	else
	{
		other_poly = NULL;
	}

	time = clock();
	batch = new Batch_smooth(lim, lbp, mfb, tree_bits, outputfile, checked_poly, param->smoothness_checking_side, param->print_limit, other_poly);
	cout << "Prime products set up, it took " << (clock() - time) / CLOCKS_PER_SEC << " seconds." << endl;

	time = clock();

	int total_num_relations = 0;

	if (inputfile.find(".bin") == string::npos)
	{
		//for reading gnfs-lasieve output.
		total_num_relations = read_batch_from_file(inputfile, batch, param);
	}
	else
	{
		//for reading binary formated files, assumed to have a .bin extension.
		total_num_relations = read_batch_from_binary_file(inputfile, batch, param);
	}

	if (total_num_relations >= 0)
	{
		cout << "Batch-checking all numbers took " << (clock() - time) / CLOCKS_PER_SEC << " seconds. Found a total of " << total_num_relations << " relations." << endl;
	}

	delete batch;
	delete param->algebraic_polynomial;
	delete param->rational_polynomial;
	delete param;

	if (total_num_relations < 0)
	{
		//The relation file could not be read.
		return -1;
	}

	cout << "Batch checking complete." << endl;
	return 0;
}
