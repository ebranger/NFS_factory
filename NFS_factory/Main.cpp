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
the batch code in CADO-NFS, and relation factoring code found in
gnfs-lasieve by Franke and Kleinjung.
----------------------------------------------------------------------*/


#include <iostream>
#include <fstream>
#include "mpirxx.h"
#include <math.h>
#include <string.h>
#include "Polynomial.h"
#include "Batch_smooth.h"
#include <time.h>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <direct.h>

using namespace std;

struct nfs_parameters {
	//only those possibly needed. In principle only one side is needed...
	Polynomial* rational_polynomial, *algebraic_polynomial;
	int lbpr, lbpa;
	int smoothness_checking_side; // 0 for rational, 1 for algebraic.
	long long maxa, maxb;
	int smooth_factor_reduce;
};

//Function for reading a polynomial. Read only the parameters of interest to this program and store in an nfs_parameter struct.
void read_polynomial(string filename, nfs_parameters* params)
{
	ifstream read;
	read.open(filename, ifstream::in);

	if (read.fail())
	{
		cout << "Failed to open input polynomial file." << endl;
		return;
	}

	//values indicating nothing read.
	params->lbpa = 0;
	params->lbpr = 0;
	params->smoothness_checking_side = 2;

	string line;
	int position;
	string data;

	mpz_class* number = new mpz_class;

	while(getline(read,line))
	{
		if (line.length() > 0 && line.substr(0,1) != "#")
		{
			position = line.find_first_of(":");
			data = line.substr(0, position);
			line = line.erase(0, position + 1);

			//For rational root given as x - m in the input. This is a short form occasionally used.
			//Otherwise, both Y0 and Y1 needs to be set. 
			if (data.length() == 1 && data[0] == 'm')
			{
				
				*number = line;
				*number = *number * -1;
				params->rational_polynomial->set_coeff(0, number);
				params->rational_polynomial->set_coeff(1, 1);
			}
			//Polynomials in GGNFS format.
			else if (data.length() == 2)
			{
				if (data[0] == 'c' && isdigit(data[1]))
				{
					*number = line;
					data = data.erase(0, 1);
					params->algebraic_polynomial->set_coeff(stoi(data), number);
				}
				else if (data[0] == 'Y' && isdigit(data[1]))
				{
					*number = line;
					data = data.erase(0, 1);
					params->rational_polynomial->set_coeff(stoi(data), number);
				}
			}
			else if (data.length() == 4)
			{
				//Maximum bits for large primes on the algebraic side
				if (data == "lpba")
				{
					params->lbpa = stoi(line);
				}
				//Maximum bits for large primes on the rational side
				else if (data == "lpbr")
				{
					params->lbpr = stoi(line);
				}
				//The side to perform batch smoothnes checking on. 0 is rational, 1 is algebraic
				else if (data == "side")
				{
					params->smoothness_checking_side = stoi(line);
				}
				//maximum a value in a relation. Relations with larger (absolute) values are discarded
				else if (data == "maxa")
				{
					params->maxa = stoi(line);
				}
				//Maximum b value in a relation. Relations with larger values are discarded
				else if (data == "maxb")
				{
					params->maxb = stoi(line);
				}
			}
			//For controlling the 2lp cutoff.
			else if (data.length() == 6)
			{
				if (data == "reduce")
				{
					params->smooth_factor_reduce = stoi(line);
				}
			}
		}
	}

	read.close();
	delete number;
}

//Function for reading relations to be processed.
int read_batch_from_file(string filename, Batch_smooth* batch)
{
	mpz_class* value = new mpz_class;

	FILE* infile;

	infile = fopen(filename.c_str(), "r");

	if (infile == NULL)
	{
		cout << "Could not open input file." << endl;
		return -1;
	}

	char buffer[512];
	long long a, b;
	bool a_negative;
	char *tmp;
	string line;
	unsigned int relation_list[20];
	unsigned int copy_testing[20];

	char *line_end;

	int relation_num = 0;

	while (!feof(infile)) {
		//read a,b,rational factor list from each line of the input.
		//fgets and char checking is so much faster than stringstream...
		fgets(buffer, sizeof(buffer), infile);

		tmp = buffer;

		if (buffer[0] == '-')
		{
			a_negative = true;
			tmp++;
		}
		else
		{
			a_negative = false;
		}

		a = 0;
		b = 0;

		//a is stored in decimal and can be negative
		while (isdigit(*tmp))
		{
			a = 10 * a + (*tmp - '0');
			tmp++;
		}

		if (a_negative)
		{
			a = a *-1;
		}

		if (*tmp != ',')
		{
			goto next_relation;
		}

		tmp++;

		//b is stored in decimal and is positive
		while (isdigit(*tmp))
		{
			b = 10 * b + (*tmp - '0');
			tmp++;
		}

		if (*tmp != ':')
		{
			goto next_relation;
		}

		*tmp++;

		if (*tmp == ':')
		{
			//empty rational list, but perhaps a complete algebraic one
			*tmp++;
		}
		line_end = strstr(tmp, ":");

		if (line_end != NULL)
		{
			//rational list, the line ends ends with ":\n\0". remove the newline.
			*line_end = '\0';
		}
		else
		{
			//Store the relation list in a char array tmp.
			line_end = strstr(tmp, "\n");

			if (line_end != NULL)
			{
				//algebraic list, the line ends ends with "\n\0". remove the newline.
				*line_end = '\0';
			}
		}//else we hope that the null termination comes at the correct position

		//Add the relation to be bacth smoothnes checked.
		batch->Add_number(a, b, tmp);
		relation_num++;

	next_relation:;
	}

	//Finish any remaining relations, even if the tree is not filled.
	batch->Do_batch_check();

	delete value;
	return batch->Get_num_relations_found();
}

//Read a list of binary values and convert to a char array for printing.
void get_relation_factor_list(unsigned int* list, int entries, char* relation_line)
{
	char* text_ptr = relation_line;

	int char_written;

	for (int i = 1; i <= list[0]; i++)
	{
		char_written = sprintf(text_ptr, "%X", list[i]);

		text_ptr += char_written;
		if (i < entries)
		{
			*text_ptr = ',';
			text_ptr++;
		}
	}

	*text_ptr = '\0';
}

//Read relations from a binary format file. 
int read_batch_from_binary_file(string filename, Batch_smooth* batch)
{
	//The binary format is (int)a, (uint)b, (int)num_primes, (int)prime1 ... (int)prime_num_primes
	char buffer[512];
	int a, b;
	bool a_negative;
	char *tmp;
	string line;
	unsigned int relation_list[20];
	int factors_in_list;
	char *line_end;

	int relation_num = 0;

	ifstream infile(filename, ios::binary);
	char* relation_line = new char[300];

	while (infile.peek() != EOF)
	{
		//read (a,b) value of relation
		infile.read(reinterpret_cast<char*>(&a), sizeof(int));
		infile.read(reinterpret_cast<char*>(&b), sizeof(unsigned int));

		//read how many primes are stored in the relation.
		infile.read(reinterpret_cast<char*>(&relation_list[0]), sizeof(int));

		//read the primes.
		for (int i = 1; i <= relation_list[0]; i++)
		{
			infile.read(reinterpret_cast<char*>(&relation_list[i]), sizeof(int));
		}
		get_relation_factor_list(relation_list, relation_list[0], relation_line);

		//Add the relation to 
		batch->Add_number(a, b, relation_line);
		relation_num++;

	}
	infile.close();

	//Check final relations, even if the tree is not filled.
	batch->Do_batch_check();

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
	param->maxa = 2147483647; //31 bit as default for now. The data I have should be filtered to have this.
	param->maxb = 2147483647;

	param->algebraic_polynomial = new Polynomial();
	param->rational_polynomial = new Polynomial();

	param->smooth_factor_reduce = 4; //i.e. make prime product up to 2^lbp{a/r} / smooth_factor_reduce, and check the remaining part up to 2^lbp{a/r} using a 2lp approach. Defualt 4, lose ~ 3% relations but reduce memory by a factor of ~3.

	read_polynomial(polyfile, param);

	//currently, lbpr/a of 31 is the maximum that can be handled. Several optimizations are needed before 32 or larger will work well.
	if (param->lbpr > 31 && param->smoothness_checking_side == 0)
	{
		cout << "Rational primes larger than 31 bits currently not supported" << endl;
		return -1;
	}

	if (param->lbpa > 31 && param->smoothness_checking_side == 1)
	{
		cout << "Algebraic primes larger than 31 bits currently not supported" << endl;
		return -1;
	}
		

	int tree_bits;

	tree_bits = 30 * param->algebraic_polynomial->get_degree() + 30; // Assume maxumum a,b is 30 bits (should be filtered to <29 bits in my SNFS-220 data....) and maximum coefficient is 30 bits (at bigger coefficients, batch factoring becomes too slow, and regular SNFS may be preferred? )

	//Check that we read the input polynomial correctly.
	if (param->algebraic_polynomial->get_degree() > 0 && param->rational_polynomial->get_degree() > 0 && param->lbpa > 0 && param->lbpr > 0 && param->smoothness_checking_side != 2)
	{
		time = clock();
		if (param->smoothness_checking_side == 0)
		{
			batch = new Batch_smooth(param->lbpr, tree_bits, outputfile, param->rational_polynomial,0, param->smooth_factor_reduce); 
		}
		else if (param->smoothness_checking_side == 1)
		{
			batch = new Batch_smooth(param->lbpa, tree_bits, outputfile, param->algebraic_polynomial,1, param->smooth_factor_reduce); 
		}
		else
		{
			cout << "Could not find which side to check for smoothness. Check the input polynomial." << endl; 
			return -1;
		}

		cout << "Prime products set up, it took " << (clock() - time) / 1000 << " seconds." << endl;
	}
	else
	{
		cout << "Failed to read input polynomial file." << endl;
		return -1;
	}

	time = clock();

	int total_num_relations = 0;

	if (inputfile.find(".bin") == string::npos)
	{
		//for reading gnfs-lasieve output.
		total_num_relations = read_batch_from_file(inputfile, batch);
	}
	else
	{
		//for reading binary formated files, assumed to have a .bin extension.
		total_num_relations = read_batch_from_binary_file(inputfile, batch);
	}

	cout << "Batch-checking all numbers took " << (clock() - time) / 1000 << " seconds. Found a total of " << total_num_relations << " relations." << endl;

	cout << "Batch checking complete." << endl;

	return 0;
}