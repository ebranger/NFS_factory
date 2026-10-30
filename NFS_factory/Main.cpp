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
and the batch factorization code in CADO-NFS.
----------------------------------------------------------------------*/


#include <iostream>
#include <fstream>
#include "mpirxx.h"
#include <math.h>
#include <string.h>
#include "Polynomial.h"
#include "Batch_smooth.h"
#include "Profile.h"
#include <time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>

using namespace std;

#define MAX_LINE_LENGTH 1024 //Longer lines in a relation file are skipped. gnfs-lasieve relations are typically less than 200 characters.
#define MAX_BINARY_PRIMES 32 //Maximum number of primes in one relation in a binary relation file. More than this is treated as a corrupt file.
#define READ_BUFFER_SIZE (16 * 1024 * 1024) //Binary relation files are read in blocks of this size.
#define STAGE2_ONLY_TREE_LEVELS 22 //Default maximum number of levels in the remainder tree of a run on a 3LP candidate file. Such runs are few, so they can use more memory than the runs on relation files.

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
	int tree_levels; // Maximum number of levels in the remainder tree. 0 means the default.
	int tree_levels_3lp; // Maximum number of levels in the remainder tree of a run on a 3LP candidate file.
	bool save_3lp; // 3LP candidates (mfb more than 2 lpb) are saved to a file instead of checked in a second stage during the run.
	string file_3lp; // File for the saved 3LP candidates. Empty means the default name.
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
			//Maximum bits for the remainder on the algebraic side after batch smoothness checking, to be factored by rho/QS
			else if (data == "mfba")
			{
				params->mfba = stoi(line);
			}
			//Maximum bits for the remainder on the rational side after batch smoothness checking, to be factored by rho/QS
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
			//Maximum number of levels in the remainder tree. Each level less halves the memory use of the tree.
			else if (data == "treelevels")
			{
				params->tree_levels = stoi(line);
			}
			//The same for runs on a 3LP candidate file.
			else if (data == "3lptreelevels")
			{
				params->tree_levels_3lp = stoi(line);
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
			//3LP candidates: checked in a second stage during the run (auto, default) or saved to a file (save).
			else if (data == "3lp")
			{
				if (line == "save")
				{
					params->save_3lp = true;
				}
				else if (line == "auto")
				{
					params->save_3lp = false;
				}
				else
				{
					cout << "Unknown 3lp setting " << line << " on line " << line_number << " of the polynomial file. Use auto or save." << endl;
					return false;
				}
			}
			//File for saved 3LP candidates.
			else if (data == "3lpfile")
			{
				params->file_3lp = line;
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
long long read_batch_from_file(string filename, Batch_smooth* batch, nfs_parameters* param)
{
	FILE* infile;

	infile = fopen(filename.c_str(), "r");

	if (infile == NULL)
	{
		cout << "Could not open input file." << endl;
		return -1;
	}
	setvbuf(infile, NULL, _IOFBF, 1 << 20);

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
			while ((c = _fgetc_nolock(infile)) != EOF && c != '\n');
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
	batch->Finish();

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

//A buffer for reading a binary file in large blocks.
struct read_buffer
{
	FILE* file;
	vector<unsigned char> data;
	size_t position; //first unread byte
	size_t filled; //bytes of data in the buffer

	//Make sure at least bytes unread bytes are in the buffer. Returns false at the end of the file.
	bool ensure(size_t bytes)
	{
		while (filled - position < bytes)
		{
			//Move the unread part to the start and fill up the rest.
			memmove(data.data(), data.data() + position, filled - position);
			filled -= position;
			position = 0;

			size_t read = fread(data.data() + filled, 1, data.size() - filled, file);
			if (read == 0)
			{
				return false;
			}
			filled += read;
		}
		return true;
	}
};

//Read relations from a binary format file. With candidates, the file holds 3LP candidates, for a batch that is the 3LP second stage.
long long read_batch_from_binary_file(string filename, Batch_smooth* batch, nfs_parameters* param, bool candidates)
{
	//The binary format is (int)a, (uint)b, (uint)num_primes, (uint)prime1 ... (uint)prime_num_primes
	//3LP candidates are followed by (uint64)low and (uint64)high 64 bits of the cofactor.
	int a;
	unsigned int b;
	unsigned int num_primes;
	unsigned int primes[MAX_SAVED_PRIMES];
	unsigned int max_primes = candidates ? MAX_SAVED_PRIMES : MAX_BINARY_PRIMES;
	size_t extra_bytes = candidates ? 2 * sizeof(unsigned long long) : 0;
	unsigned long long cofactor_limbs[2];
	mpz_class cofactor;

	long long relation_num = 0;
	long long skipped_relations = 0;
	bool corrupt = false;
	bool stopped = false;

	read_buffer input;
	input.file = fopen(filename.c_str(), "rb");

	if (input.file == NULL)
	{
		cout << "Failed to open input saved relations file." << endl;
		return -1;
	}

	input.data.resize(READ_BUFFER_SIZE);
	input.position = 0;
	input.filled = 0;

	while (input.ensure(3 * sizeof(unsigned int)))
	{
		//read (a,b) value of relation, and how many primes are stored in the relation.
		const unsigned char* record = input.data.data() + input.position;
		memcpy(&a, record, sizeof(int));
		memcpy(&b, record + 4, sizeof(unsigned int));
		memcpy(&num_primes, record + 8, sizeof(unsigned int));

		if (num_primes > max_primes || !input.ensure((3 + num_primes) * sizeof(unsigned int) + extra_bytes))
		{
			corrupt = true;
			break;
		}

		//read the primes.
		record = input.data.data() + input.position;
		memcpy(primes, record + 12, num_primes * sizeof(unsigned int));
		if (candidates)
		{
			memcpy(cofactor_limbs, record + (3 + num_primes) * sizeof(unsigned int), sizeof(cofactor_limbs));
		}
		input.position += (3 + num_primes) * sizeof(unsigned int) + extra_bytes;

		if (!relation_in_range(a, b, param))
		{
			skipped_relations++;
			continue;
		}

		//Add the relation to be batch smoothness checked.
		if (candidates)
		{
			mpz_set_ui(cofactor.get_mpz_t(), (mpir_ui)cofactor_limbs[1]);
			mpz_mul_2exp(cofactor.get_mpz_t(), cofactor.get_mpz_t(), 64);
			mpz_add_ui(cofactor.get_mpz_t(), cofactor.get_mpz_t(), (mpir_ui)cofactor_limbs[0]);
			batch->Add_candidate(a, b, primes, num_primes, cofactor.get_mpz_t());
		}
		else
		{
			batch->Add_number_binary(a, b, primes, num_primes);
		}
		relation_num++;

		if (param->maxrels > 0 && relation_num >= param->maxrels)
		{
			cout << "Stopping after " << relation_num << " relations (maxrels)." << endl;
			stopped = true;
			break;
		}
	}

	if (!corrupt && !stopped && input.filled > input.position)
	{
		//Part of a relation is left at the end of the file.
		corrupt = true;
	}
	fclose(input.file);

	//Check final relations, even if the tree is not filled.
	batch->Finish();

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
	unsigned long long program_start = profile_now();

	nfs_parameters* param = new nfs_parameters();
	param->maxa = LLONG_MAX; //No limit unless set in the polynomial file.
	param->maxb = LLONG_MAX;
	param->maxrels = 0;

	param->print_limit = 1000; //Msieve default.
	param->cado_format = false;
	param->tree_levels = 0;
	param->tree_levels_3lp = STAGE2_ONLY_TREE_LEVELS;
	param->save_3lp = false;

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

	//Cofactors of more than 2 lbp bits can have three large primes, more than that is not supported.
	if (lbp * 3 < mfb)
	{
		cout << "More than three large primes is not supported. Use mfb" << side_name << " of at most 3 lpb" << side_name << "." << endl;
		return -1;
	}

	//A file of 3LP candidates is checked in the 3LP second stage only.
	bool candidate_input = inputfile.size() > 4 && inputfile.compare(inputfile.size() - 4, 4, ".3lp") == 0;
	bool three_lp = mfb > 2 * lbp && !candidate_input;

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
	unsigned long long t = profile_now();
	if (candidate_input)
	{
		//The cofactors are less than 2^(3 lbp). The tree size is limited by 3lptreelevels instead of treelevels.
		cout << "Checking 3LP candidates against all primes up to 2^" << lbp << "." << endl;
		batch = new Batch_smooth(lim, lbp, mfb, 3 * lbp, outputfile, checked_poly, param->smoothness_checking_side, param->print_limit, other_poly, param->tree_levels_3lp, true);
	}
	else
	{
		batch = new Batch_smooth(lim, lbp, mfb, tree_bits, outputfile, checked_poly, param->smoothness_checking_side, param->print_limit, other_poly, param->tree_levels);
	}

	if (three_lp)
	{
		Batch_smooth* second_stage = NULL;
		string candidate_filename;

		if (param->save_3lp)
		{
			//By default, one candidate file for each output file and input file, so that runs on different input files can run at the same time.
			candidate_filename = param->file_3lp;
			if (candidate_filename.empty())
			{
				size_t name_start = inputfile.find_last_of("\\/");
				candidate_filename = outputfile + "." + inputfile.substr((name_start == string::npos) ? 0 : name_start + 1) + ".3lp";
			}
			cout << "Saving 3LP candidates of up to " << mfb << " bits to " << candidate_filename << endl;
		}
		else
		{
			cout << "Checking 3LP candidates of up to " << mfb << " bits against all primes up to 2^" << lbp << " during the run." << endl;
			second_stage = new Batch_smooth(lim, lbp, mfb, 3 * lbp, outputfile, checked_poly, param->smoothness_checking_side, param->print_limit, other_poly, param->tree_levels, true);
		}

		if (!batch->Set_three_lp_handling(second_stage, candidate_filename))
		{
			cout << "Could not open the 3LP candidate file " << candidate_filename << endl;
			delete batch;
			return -1;
		}
	}
	profile_add(PROF_SETUP, t);
	cout << "Prime products set up, it took " << (clock() - time) / CLOCKS_PER_SEC << " seconds." << endl;

	time = clock();

	long long total_num_relations = 0;

	//The time in the readers that is not spent on batch checking is counted as reading.
	unsigned long long batch_cycles = 0;
	for (int i = PROF_POLY_EVAL; i < PROF_NUM_PHASES; i++)
	{
		batch_cycles -= profile_cycles[i];
	}
	t = profile_now();

	if (candidate_input)
	{
		//3LP candidates saved by an earlier run, with a .3lp extension.
		total_num_relations = read_batch_from_binary_file(inputfile, batch, param, true);
	}
	else if (inputfile.find(".bin") == string::npos)
	{
		//for reading gnfs-lasieve output.
		total_num_relations = read_batch_from_file(inputfile, batch, param);
	}
	else
	{
		//for reading binary formated files, assumed to have a .bin extension.
		total_num_relations = read_batch_from_binary_file(inputfile, batch, param, false);
	}

	unsigned long long reader_cycles = profile_now() - t;
	for (int i = PROF_POLY_EVAL; i < PROF_NUM_PHASES; i++)
	{
		batch_cycles += profile_cycles[i];
	}
	profile_cycles[PROF_READ] += reader_cycles - batch_cycles;

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

	print_profile(profile_now() - program_start);
	cout << "Batch checking complete." << endl;
	return 0;
}
