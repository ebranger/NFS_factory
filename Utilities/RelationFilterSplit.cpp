// A small program to split a relation file into smaller files based on the largest prime in the relation.

//TODO: allow setting startvalue, endvalue, number_of_bins at program startup rather than hard-coding them.

#include "mpirxx.h"
#include "Polynomial.h"

#include <fstream>
#include <iostream>
#include <string>
#include <cerrno>
#include <sstream>
#include <algorithm>
#include <math.h>
using namespace std;


double startvalue = 10000000.0; // minimum spq
double endvalue = 536870912.0; // with lpbr = 29
const int number_of_bins = 360;
long max_value_of_a_b = 536870912; //2^29

long valid_relation, discarded_relations, large_relations, bad_relations;

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

	while (getline(read, line))
	{
		if (line.length() > 0 && line.substr(0, 1) != "#")
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

int convert_relation_text(char* text, unsigned int* list)
{
	list[0] = 0;

	unsigned int digit = 0;

	while (*text != ':' && *text != '\n' && *text != '\0')
	{

		if (*text >= '0' && *text <= '9')
		{
			digit = digit * 16 + *text - '0';
			text++;
		}
		else if (*text >= 'a' && *text <= 'f')
		{
			digit = digit * 16 + *text - 'a' + 10;
			text++;
		}
		else if (*text >= 'A' && *text <= 'F')
		{
			digit = digit * 16 + *text - 'A' + 10;
			text++;
		}
		else if (*text == ',' || *text == ':' || *text == '\n')
		{
			list[0]++;
			list[list[0]] = digit;
			digit = 0;
			text++;
		}
		else
		{
			return -1; // could not parse relation
		}

	}
	if (digit != 0)
	{
		list[0]++;
		list[list[0]] = digit;
	}


	return 0;
}

void get_relation_line(long long a, long long b, unsigned int* list, int entries, char* relation_line, string side_to_test)
{
	char* text_ptr = relation_line;

	int char_written;

	char_written = sprintf(text_ptr, "%lld", a);

	text_ptr += char_written;
	*text_ptr = ',';
	text_ptr++;

	char_written = sprintf(text_ptr, "%lld", b);

	text_ptr += char_written;
	*text_ptr = ':';
	text_ptr++;

	if (side_to_test == "a")
	{
		*text_ptr = ':';
		text_ptr++;
	}
	*text_ptr = '\0';

	for (int i = 1; i <= list[0]; i++)
	{
		char_written = sprintf(text_ptr, "%x", list[i]);

		text_ptr += char_written;
		if (i < entries)
		{
			*text_ptr = ',';
			text_ptr++;
		}
	}

	if (side_to_test == "r")
	{
		*text_ptr = ':';
		text_ptr++;
	}
	*text_ptr = '\n';
	text_ptr++;
	*text_ptr = '\0';
}

int get_bin(int number)
{
	//Bin division is given by global constants defined at the beginning of this file.

	if (number == 0)
	{
		//Will happen if ther are no printed primes in the relation, i.e. all primes are below 1000. Just put those in the first bin.
		return 1;
	}

	double temp;

	//Change: use a linear bining instead. 

	temp = floor(((number - startvalue) / (endvalue - startvalue))*number_of_bins); // should be in the range 1 to 50.

	temp = temp + 1;

	if (temp <= 0)
	{
		cout << "Found a relation with the maximum prime below " << startvalue << ", but sieving was above that value..." << endl;
		temp = 1;
	}
	if (temp > number_of_bins)
	{
		cout << "Found a relation with the maximum prime above " << endvalue << ", but a lower maximmum value was specified..." << endl;
		temp = number_of_bins;
	}

	return temp;
}

bool is_valid_relation(long long a, long long b, unsigned int* list, int entries, Polynomial* poly, mpz_class* temp_bigint)
{
	if (abs(a) > max_value_of_a_b || b > max_value_of_a_b)
	{
		large_relations++;
		return false;
	}

	poly->evaluate_polynomial_homogenous(temp_bigint, a, b);

	for (int i = 1; i <= entries; i++) // note: list[0] = entries, followed by the list of primes starting at list[1].
	{
		if (*temp_bigint % list[i] != 0)
		{
			bad_relations++;
			return false;
		}
		else
		{
			*temp_bigint = *temp_bigint / list[i];
		}
	}

	return true;

}

int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		cout << "Usage: RelationFileSplit <input file name> <output file name without .out> <polynomial file> <Side that was sieved, r or a>" << endl;
		return 0;
	}

	string inputfile, outputfile, polyfile, side_to_test, currentoutput;
	inputfile = argv[1];
	outputfile = argv[2];
	polyfile = argv[3];
	side_to_test = argv[4];

	valid_relation = 0;
	discarded_relations = 0;
	large_relations = 0;
	bad_relations = 0;


	Polynomial* poly;
	mpz_class* temp_bigint = new mpz_class;

	nfs_parameters* param = new nfs_parameters();
	param->algebraic_polynomial = new Polynomial();
	param->rational_polynomial = new Polynomial();
	read_polynomial(polyfile, param);

	bool line_is_ok;

	if (side_to_test == "r")
	{
		poly = param->rational_polynomial;
	}
	else if (side_to_test == "a")
	{
		poly = param->algebraic_polynomial;
	}
	else
	{
		cout << "The side that was sieved must be specified as r or a." << endl;
		return -1;
	}

	//string relations[65];
	int batch_write_relations[number_of_bins + 2];

	char** relations = new char*[number_of_bins + 2];
	char** write_here = new char*[number_of_bins + 2];

	for (int i = 0; i < number_of_bins + 2; i++)
	{
		relations[i] = new char[1500000];
		batch_write_relations[i] = 0;
		write_here[i] = relations[i];
	}

	stringstream tohex;
	tohex << hex;
	//string relation_line;
	char* relation_line = new char[150];

	FILE* infile, *outfile;

	infile = fopen(inputfile.c_str(), "r");

	if (infile == NULL)
	{
		cout << "Could not open relation input file." << endl;
		return -1;
	}

	char buffer[512];
	long long a, b;
	bool a_negative;
	char *tmp;
	string line;
	unsigned int relation_list[20];
	int factors_in_list;
	int bin = 0;

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

		while (isdigit(*tmp))
		{
			a = 10 * a + (*tmp - '0');
			tmp++;
		}

		if (a_negative)
		{
			a = a * -1;
		}

		if (*tmp != ',')
		{
			goto next_relation;
		}

		tmp++;

		while (isdigit(*tmp))
		{
			b = 10 * b + (*tmp - '0');
			tmp++;
		}

		if (*tmp != ':')
		{
			goto next_relation;
		}

		tmp++;
		if (side_to_test == "r")
		{
			if (*tmp == ':')
			{
				//empty rational list, but perhaps a complete algebraic one
				*tmp = '\0';
			}
			else
			{
				line_end = strstr(tmp, ":");
				if (line_end != NULL)
				{
					//rational list, the line ends ends with ":\n\0". remove the newline.
					*line_end = '\0';
				}
				else
				{
					//Malformed relation?
					goto next_relation;
				}
			}
		}
		else if (side_to_test == "a")
		{
			tmp = strstr(tmp, ":");
			tmp++;

			line_end = strstr(tmp, "\n");
			if (line_end != NULL)
			{
				//rational list, the line ends ends with ":\n\0". remove the newline.
				*line_end = '\0';
			}
			else
			{
				//Malformed relation?
				goto next_relation;
			}
		}


		convert_relation_text(tmp, relation_list);

		factors_in_list = relation_list[0];

		sort(relation_list + 1, relation_list + 1 + factors_in_list);

		get_relation_line(a, b, relation_list, factors_in_list, relation_line, side_to_test);

		if (is_valid_relation(a, b, relation_list, factors_in_list, poly, temp_bigint))
		{
			bin = get_bin(relation_list[factors_in_list]);

			strcpy(write_here[bin], relation_line);

			write_here[bin] += strlen(relation_line);

			//relations[bin] += relation_line;
			batch_write_relations[bin]++;

			valid_relation++;
		}
		else
		{
			discarded_relations++;
		}



		if (batch_write_relations[bin] > 10000)
		{
			currentoutput = outputfile + to_string(bin) + ".out";
			outfile = fopen(currentoutput.c_str(), "a");

			if (outfile == NULL)
			{
				cout << "Could not open output file." << endl;
				return -1;
			}
			fputs(relations[bin], outfile);
			fclose(outfile);
			batch_write_relations[bin] = 0;
			//relations[bin] = "";
			write_here[bin] = relations[bin];
			*write_here[bin] = '\0';
		}

		//Add to list, save list when too long

		relation_num++;
		if (relation_num % 1000000 == 0)
		{
			cout << "Filtering has handled " << relation_num / 1000000 << " M relations." << endl;
		}

	next_relation:;
	}

	for (int i = 1; i <= number_of_bins; i++)
	{
		if (batch_write_relations[i] > 0)
		{
			currentoutput = outputfile + to_string(i) + ".out";
			outfile = fopen(currentoutput.c_str(), "a");

			if (outfile == NULL)
			{
				cout << "Could not open output file." << endl;
				return -1;
			}
			fputs(relations[i], outfile);
			fclose(outfile);
		}
	}

	delete poly;
	delete temp_bigint;

	cout << "Found " << valid_relation << " valid relations, and discarded " << discarded_relations << " relations. ( " << large_relations << " large relations, " << bad_relations << " bad relations.)" << endl;

	return 0;
}