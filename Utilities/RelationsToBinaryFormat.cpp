// A program to convert a relation foel from GGNFS uncompressed output to a binary format or the other way around.
//
// TODO: allow the program to also convert relations that are for the algebraic side, not just the rational one. 

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

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

	list[0]++;
	list[list[0]] = digit;

	return 0;
}

void get_relation_text_line(int a, int b, unsigned int* list, int entries, char* relation_line)
{
	char* text_ptr = relation_line;

	int char_written;

	char_written = sprintf(text_ptr, "%d", a);

	text_ptr += char_written;
	*text_ptr = ',';
	text_ptr++;

	char_written = sprintf(text_ptr, "%d", b);

	text_ptr += char_written;
	*text_ptr = ':';
	text_ptr++;
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

	*text_ptr = ':';
	text_ptr++;
	*text_ptr = '\n';
	text_ptr++;
	*text_ptr = '\0';
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		cout << "Usage: ReltationsToBinaryFormat <input file name> <output file name> <b/t (output to binary or text)>" << endl;
		return 0;
	}

	string inputfile, outputfile, mode;
	inputfile = argv[1];
	outputfile = argv[2];
	mode = argv[3];

	char buffer[512];
	int a, b;
	bool a_negative;
	char *tmp;
	string line;
	unsigned int relation_list[20];
	int factors_in_list;
	char *line_end;

	int num_relations = 0;

	if (mode[0] == 'b')
	{

		FILE* infile, *outfile;

		infile = fopen(inputfile.c_str(), "r");

		if (infile == NULL)
		{
			cout << "Could not open input file." << endl;
			return -1;
		}

		outfile = fopen(outputfile.c_str(), "wb");

		if (outfile == NULL)
		{
			cout << "Could not open output file." << endl;
			return -1;
		}

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
				line_end = strstr(tmp, "\n");

				if (line_end != NULL)
				{
					//algebraic list, the line ends ends with "\n\0". remove the newline.
					*line_end = '\0';
				}



			}//else we hope that the null termination comes at the correct position

			convert_relation_text(tmp, relation_list);

			fwrite(&a, sizeof(unsigned int), 1, outfile);
			fwrite(&b, sizeof(int), 1, outfile);
			fwrite(relation_list, sizeof(unsigned int), relation_list[0] + 1, outfile);
			num_relations++;
		next_relation:;
			if (num_relations % 1000000 == 0)
			{
				cout << "Converted " << num_relations << " relations." << endl;
			}
		}
		fclose(outfile);
		fclose(infile);
	}
	else if (mode[0] == 't')
	{
		ifstream infile(inputfile, ios::binary);
		ofstream outfile(outputfile);
		char* relation_line = new char[300];

		while (infile.peek() != EOF)
		{
			infile.read(reinterpret_cast<char*>(&a), sizeof(int));
			infile.read(reinterpret_cast<char*>(&b), sizeof(unsigned int));
			infile.read(reinterpret_cast<char*>(&relation_list[0]), sizeof(int));

			for (int i = 1; i <= relation_list[0]; i++)
			{
				infile.read(reinterpret_cast<char*>(&relation_list[i]), sizeof(int));
			}
			get_relation_text_line(a, b, relation_list, relation_list[0], relation_line);
			outfile << relation_line;
			num_relations++;
			if (num_relations % 1000000 == 0)
			{
				cout << "Converted " << num_relations << " relations." << endl;
			}
		}
		infile.close();
		outfile.close();
	}
	else
	{
		cout << "unknown conversion mode: " << mode << ", should be b or t." << endl;
	}

	

	cout << "Converted " << num_relations << " relations." << endl;
	return 0;
}

