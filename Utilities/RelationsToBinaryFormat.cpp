// A program to convert a relation file from GGNFS uncompressed output to a binary format, the other way around, or to a binary
// format of only a and b (.ab) for NFS_factory.
//
// Usage: RelationsToBinaryFormat <input file> <output file> b [side]
//        RelationsToBinaryFormat <input file> <output file> a [side]
//        RelationsToBinaryFormat <input file> <output file> t <side> <polynomial file>
//
// Modes:
//   b  text to binary (.bin): (int)a, (uint)b, (uint)number of primes, (uint)primes of the side that was sieved.
//   t  binary (.bin) to text, with the primes in the list of the side that was sieved (the list of the other side is empty).
//   a  text, or binary if the input file name ends with .bin, to only a and b (.ab): a header of 16 bytes ("NFAB", (uint)version 1,
//      (uint)side whose primes the relations were saved for, (uint)flags 0), then (int64)a, (uint32)b per relation. NFS_factory
//      rebuilds the factor list of that side for the relations it finds.
//
// side is the side that the relations were sieved on, whose factor list is kept: 0 or r for rational (the default), 1 or a for
// algebraic. The .bin format does not say which side its primes are for (NFS_factory takes them for the side it does not batch
// check), so the same side must be given when a .bin file is converted on.
//
// Before converting, the first relations of the input are checked against the side, and a file that fails the check is not
// converted. In a text file, the lines with only one factor list (as RelationFilterSplit writes them) must mostly have it on that
// side; a relation can have an empty list, so a few lines are looked at. Text with both lists, as the siever writes it, has the
// primes of both sides, and either side can be kept. For t, the primes of the .bin file must divide the norms on that side, from the
// polynomial file, since the .bin file does not record its side. A .bin file converted to .ab with the wrong side is refused by
// NFS_factory, which checks the side in the header of the .ab file against the side in its polynomial file.
//
// The binary format holds a in 32 bits, so relations with |a| of 2^31 or more cannot be converted to it, and are skipped. Binary files
// written by older versions of this program can hold such relations with a wrong a (a changed by 2^32): those do not give valid
// relations, and NFS_factory leaves them out.

#define _CRT_SECURE_NO_WARNINGS
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

static const unsigned int MAX_PRIMES = 64; //primes of one side of a relation
static const size_t LINE_SIZE = 1 << 16; //longest text line read
static const int CHECK_RELATIONS = 100; //relations checked against the side before converting

struct Relation
{
	long long a, b;
	unsigned int primes[MAX_PRIMES]; //of the side that was sieved
	unsigned int num_primes;
	unsigned int other_primes; //how many the list of the other side has (text input), or 0
};

// Parse a decimal integer with an optional minus sign. Returns a pointer to the character after it, or NULL if there is none or it
// has more than 18 digits.
static const char* parse_integer(const char* text, long long* value)
{
	bool negative = false;
	int digits = 0;
	if (*text == '-')
	{
		negative = true;
		text++;
	}
	*value = 0;
	while (isdigit((unsigned char)*text))
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

// Parse a factor list (hexadecimal primes below 2^32, separated by commas, possibly empty) that ends at stop (':' or the end of the
// line). The primes go to primes if it is not NULL (at most MAX_PRIMES of them). Returns a pointer to the character that ends the
// list, or NULL if the list is malformed.
static const char* parse_list(const char* text, char stop, unsigned int* primes, unsigned int* num_primes)
{
	unsigned int n = 0;
	while (*text != stop)
	{
		char* end;
		const unsigned long long p = strtoull(text, &end, 16);
		if (end == text || !isxdigit((unsigned char)*text) || p > 0xFFFFFFFFULL || p < 2 || (primes != NULL && n >= MAX_PRIMES))
		{
			return NULL;
		}
		if (primes != NULL)
		{
			primes[n] = (unsigned int)p;
		}
		n++;
		text = end;
		if (*text == ',' && text[1] != stop)
		{
			text++;
		}
		else if (*text != stop)
		{
			return NULL;
		}
	}
	*num_primes = n;
	return text;
}

// Parse a relation line "a,b:rational primes:algebraic primes" (without its line ending). The primes of the side that was sieved go
// to r. Returns false if the line is not a relation.
static bool parse_relation(const char* line, unsigned int side, Relation* r)
{
	const char* text = parse_integer(line, &r->a);
	if (text == NULL || *text != ',')
	{
		return false;
	}
	text = parse_integer(text + 1, &r->b);
	if (text == NULL || *text != ':' || r->b <= 0 || r->b >= (1LL << 32))
	{
		return false;
	}
	unsigned int rational_count, algebraic_count;
	text = parse_list(text + 1, ':', (side == 0) ? r->primes : NULL, &rational_count);
	if (text == NULL)
	{
		return false;
	}
	text = parse_list(text + 1, '\0', (side == 1) ? r->primes : NULL, &algebraic_count);
	if (text == NULL)
	{
		return false;
	}
	r->num_primes = (side == 0) ? rational_count : algebraic_count;
	r->other_primes = (side == 0) ? algebraic_count : rational_count;
	return true;
}

static bool ends_with(const string& name, const char* suffix)
{
	const size_t n = strlen(suffix);
	return name.size() >= n && name.compare(name.size() - n, n, suffix) == 0;
}

static const char* side_name(unsigned int side)
{
	return (side == 0) ? "rational" : "algebraic";
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------------------------------------------------------------

class Reader
{
public:
	enum Result { relation, not_relation, end, error };

	Reader(bool binary, unsigned int side) : binary(binary), side(side), line(LINE_SIZE) {}

	~Reader()
	{
		if (file != NULL)
		{
			fclose(file);
		}
	}

	bool Open(const string& name)
	{
		file = fopen(name.c_str(), "rb");
		if (file != NULL)
		{
			setvbuf(file, NULL, _IOFBF, 1 << 22);
		}
		return file != NULL;
	}

	bool Rewind()
	{
		incomplete_last_line = false;
		return fseek(file, 0, SEEK_SET) == 0;
	}

	Result Read(Relation* r)
	{
		return binary ? Read_binary(r) : Read_text(r);
	}

	bool incomplete_last_line = false; //a text file that ends in the middle of a line, as a siever that is stopped may leave it
	string message; //what went wrong, for error

private:
	Result Read_text(Relation* r)
	{
		while (fgets(line.data(), (int)line.size(), file) != NULL)
		{
			size_t length = strlen(line.data());
			if (length == 0 || line[length - 1] != '\n')
			{
				if (feof(file))
				{
					incomplete_last_line = true;
					break;
				}
				//Longer than any relation: skip the rest of it.
				int c;
				while ((c = getc(file)) != EOF && c != '\n')
				{
				}
				return not_relation;
			}
			line[--length] = '\0';
			if (length > 0 && line[length - 1] == '\r')
			{
				line[--length] = '\0';
			}
			if (length == 0 || line[0] == '#')
			{
				continue;
			}
			return parse_relation(line.data(), side, r) ? relation : not_relation;
		}
		if (ferror(file))
		{
			message = "Could not read the input file.";
			return error;
		}
		return end;
	}

	// (int)a, (uint)b, (uint)number of primes, (uint)primes
	Result Read_binary(Relation* r)
	{
		unsigned int head[3];
		const size_t got = fread(head, sizeof(unsigned int), 3, file);
		if (got == 0 && feof(file))
		{
			return end;
		}
		if (got == 3 && head[2] > MAX_PRIMES)
		{
			message = "The binary file is corrupt: a relation with " + to_string(head[2]) + " primes.";
			return error;
		}
		if (got != 3 || fread(r->primes, sizeof(unsigned int), head[2], file) != head[2])
		{
			message = ferror(file) ? "Could not read the input file." : "The binary file is truncated (it ends in the middle of a relation).";
			return error;
		}
		r->a = (int)head[0];
		r->b = head[1];
		r->num_primes = head[2];
		r->other_primes = 0;
		return relation;
	}

	bool binary;
	unsigned int side;
	FILE* file = NULL;
	vector<char> line;
};

// ---------------------------------------------------------------------------------------------------------------------------------
// The check of the side
// ---------------------------------------------------------------------------------------------------------------------------------

// The polynomials of a polynomial file, as decimal strings: c0, c1, ... and Y0, Y1 (or m, for x - m).
struct Polynomials
{
	vector<string> coefficients[2]; //[side][power]
};

static bool read_polynomials(const string& name, Polynomials* poly)
{
	ifstream in(name);
	if (!in)
	{
		cout << "Could not open the polynomial file " << name << "." << endl;
		return false;
	}
	string line;
	while (getline(in, line))
	{
		const size_t colon = line.find(':');
		if (line.empty() || line[0] == '#' || colon == string::npos)
		{
			continue;
		}
		string key = line.substr(0, colon), value = line.substr(colon + 1);
		key.erase(key.find_last_not_of(" \t") + 1);
		value.erase(0, value.find_first_not_of(" \t"));
		value.erase(value.find_last_not_of(" \t\r") + 1);
		int s = -1, power = 0;
		if (key == "m" && !value.empty())
		{
			poly->coefficients[0] = { (value[0] == '-') ? value.substr(1) : "-" + value, "1" };
			continue;
		}
		if (key.size() == 2 && isdigit((unsigned char)key[1]) && (key[0] == 'c' || key[0] == 'Y'))
		{
			s = (key[0] == 'c') ? 1 : 0;
			power = key[1] - '0';
		}
		if (s < 0)
		{
			continue;
		}
		const size_t start = (!value.empty() && value[0] == '-') ? 1 : 0;
		if (value.size() == start || value.find_first_not_of("0123456789", start) != string::npos)
		{
			cout << "The value of " << key << " in the polynomial file is not a number: " << value << endl;
			return false;
		}
		vector<string>& c = poly->coefficients[s];
		if ((int)c.size() <= power)
		{
			c.resize(power + 1, "0");
		}
		c[power] = value;
	}
	for (int s = 0; s < 2; s++)
	{
		if (poly->coefficients[s].size() < 2)
		{
			cout << "The polynomial file has no " << side_name(s) << " polynomial." << endl;
			return false;
		}
	}
	return true;
}

static uint64_t residue(const string& value, uint64_t p)
{
	uint64_t r = 0;
	for (size_t i = (value[0] == '-') ? 1 : 0; i < value.size(); i++)
	{
		r = (10 * r + (uint64_t)(value[i] - '0')) % p;
	}
	return (value[0] == '-' && r != 0) ? p - r : r;
}

// Whether p divides the norm of a,b on the side with the coefficients c: the sum of c[i] a^i b^(d-i).
static bool divides_norm(const vector<string>& c, long long a, long long b, uint64_t p)
{
	const uint64_t ap = (uint64_t)(((a % (long long)p) + (long long)p) % (long long)p), bp = (uint64_t)b % p;
	uint64_t sum = 0, apow = 1;
	for (size_t i = 0; i < c.size(); i++)
	{
		uint64_t term = residue(c[i], p) * apow % p;
		for (size_t j = i + 1; j < c.size(); j++)
		{
			term = term * bp % p;
		}
		sum = (sum + term) % p;
		apow = apow * ap % p;
	}
	return sum == 0;
}

// Checks the first relations of the input against the side: for text, that the lines with only one factor list have it on that side
// (a relation can have an empty list, so it takes most of them); with poly (for .bin input), that the primes divide the norms on that
// side for most of them. Returns false, with a message, if they look sieved on the other side.
static bool check_side(Reader& reader, unsigned int side, const Polynomials* poly)
{
	int relations = 0, one_list = 0, empty_where_other_not = 0, kept_lists = 0, divide_side = 0, divide_other = 0;
	Relation r;
	for (Reader::Result result; relations < CHECK_RELATIONS && (result = reader.Read(&r)) != Reader::end; )
	{
		if (result == Reader::error)
		{
			cout << reader.message << endl;
			return false;
		}
		if (result != Reader::relation)
		{
			continue;
		}
		relations++;
		if ((r.num_primes == 0) != (r.other_primes == 0))
		{
			one_list++;
			empty_where_other_not += (r.num_primes == 0);
		}
		if (r.num_primes == 0)
		{
			continue;
		}
		kept_lists++;
		if (poly != NULL)
		{
			bool on_side = true, on_other = true;
			for (unsigned int i = 0; i < r.num_primes; i++)
			{
				on_side = on_side && divides_norm(poly->coefficients[side], r.a, r.b, r.primes[i]);
				on_other = on_other && divides_norm(poly->coefficients[1 - side], r.a, r.b, r.primes[i]);
			}
			divide_side += on_side;
			divide_other += on_other;
		}
	}
	if (!reader.Rewind())
	{
		cout << "Could not go back to the start of the input file." << endl;
		return false;
	}
	const string other = "side " + to_string(1 - side) + " (" + side_name(1 - side) + ")";
	if (2 * empty_where_other_not > one_list)
	{
		cout << "Of the first " << relations << " relations, " << one_list << " have only one factor list, and in " << empty_where_other_not
			<< " of them it is the " << side_name(1 - side) << " list: they look sieved on " << other << ". Give that side, or check the input file." << endl;
		return false;
	}
	if (poly != NULL && kept_lists > 0 && 2 * divide_side < kept_lists)
	{
		cout << "The primes kept for side " << side << " (" << side_name(side) << ") divide the " << side_name(side) << " norms of only " << divide_side
			<< " of the first " << kept_lists << " relations with primes";
		if (divide_other > divide_side)
		{
			cout << ", but the " << side_name(1 - side) << " norms of " << divide_other << ": they look sieved on " << other << ". Give that side.";
		}
		else
		{
			cout << ". Check the polynomial file and the input file.";
		}
		cout << endl;
		return false;
	}
	if (poly != NULL && kept_lists > 0)
	{
		cout << "Checked the first " << kept_lists << " relations with primes: " << divide_side << " of them divide the " << side_name(side) << " norms." << endl;
	}
	return true;
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------------------------------------------------------------

static bool write_relation(FILE* out, const string& mode, unsigned int side, const Relation& r)
{
	if (mode == "b")
	{
		const unsigned int head[3] = { (unsigned int)(int)r.a, (unsigned int)r.b, r.num_primes };
		return fwrite(head, sizeof(head), 1, out) == 1 && fwrite(r.primes, sizeof(unsigned int), r.num_primes, out) == r.num_primes;
	}
	if (mode == "a")
	{
		const unsigned int b32 = (unsigned int)r.b;
		return fwrite(&r.a, sizeof(long long), 1, out) == 1 && fwrite(&b32, sizeof(unsigned int), 1, out) == 1;
	}
	char text[64 + MAX_PRIMES * 12];
	char* p = text + sprintf(text, (side == 0) ? "%lld,%lld:" : "%lld,%lld::", r.a, r.b);
	for (unsigned int i = 0; i < r.num_primes; i++)
	{
		p += sprintf(p, (i == 0) ? "%x" : ",%x", r.primes[i]);
	}
	if (side == 0)
	{
		*p++ = ':';
	}
	*p++ = '\n';
	return fwrite(text, 1, p - text, out) == (size_t)(p - text);
}

static void usage()
{
	cout << "Usage: RelationsToBinaryFormat <input file> <output file> b [side]       text to binary (.bin)" << endl;
	cout << "       RelationsToBinaryFormat <input file> <output file> a [side]       text or binary (.bin) to only a and b (.ab)" << endl;
	cout << "       RelationsToBinaryFormat <input file> <output file> t <side> <polynomial file>   binary (.bin) to text" << endl;
	cout << "  side: the side that the relations were sieved on, whose primes are kept: 0 or r for rational (default), 1 or a for algebraic." << endl;
	cout << "  The first relations are checked against the side: in text input, the lines with only one factor list must have it on that side;" << endl;
	cout << "  for t, the primes must divide the norms on that side, from the polynomial file." << endl;
}

int main(int argc, char* argv[])
{
	const bool to_text = (argc >= 4 && string(argv[3]) == "t");
	if (argc < 4 || argc > 6 || (to_text && argc != 6) || (!to_text && argc == 6))
	{
		usage();
		if (argc == 6)
		{
			cout << "The polynomial file is only given for t." << endl;
		}
		else if (to_text)
		{
			cout << "t needs the side and the polynomial file, to check that the .bin file is for that side." << endl;
		}
		return (argc == 1) ? 0 : -1;
	}
	const string inputfile = argv[1], outputfile = argv[2], mode = argv[3];
	const string side_text = (argc >= 5) ? argv[4] : "0";
	if (mode != "b" && mode != "t" && mode != "a")
	{
		cout << "Unknown conversion mode: " << mode << ", should be b, t or a." << endl;
		return -1;
	}
	if (side_text != "0" && side_text != "1" && side_text != "r" && side_text != "a")
	{
		cout << "The side must be 0 or r (rational), or 1 or a (algebraic), not " << side_text << "." << endl;
		return -1;
	}
	const unsigned int side = (side_text == "0" || side_text == "r") ? 0 : 1;
	Polynomials poly;
	if (to_text && !read_polynomials(argv[5], &poly))
	{
		return -1;
	}

	const bool binary_input = (mode == "t") || (mode == "a" && ends_with(inputfile, ".bin"));
	if (ends_with(inputfile, ".ab"))
	{
		cout << "The input is an .ab file, which has no factor lists to convert." << endl;
		return -1;
	}
	if (mode != "t" && ends_with(inputfile, ".bin") != binary_input)
	{
		cout << "A .bin file can only be converted to text (t) or to .ab (a)." << endl;
		return -1;
	}
	Reader reader(binary_input, side);
	if (!reader.Open(inputfile))
	{
		cout << "Could not open the input file " << inputfile << "." << endl;
		return -1;
	}
	if (!check_side(reader, side, to_text ? &poly : NULL))
	{
		cout << "Nothing was converted." << endl;
		return -1;
	}

	FILE* outfile = fopen(outputfile.c_str(), "wb");
	if (outfile == NULL)
	{
		cout << "Could not open the output file " << outputfile << "." << endl;
		return -1;
	}
	setvbuf(outfile, NULL, _IOFBF, 1 << 22);
	bool ok = true;
	if (mode == "a")
	{
		unsigned int header[4] = { 0, 1, side, 0 };
		memcpy(header, "NFAB", 4);
		ok = fwrite(header, sizeof(header), 1, outfile) == 1;
	}

	long long num_relations = 0, num_failed = 0, num_too_large = 0;
	Relation r;
	for (Reader::Result result; ok && (result = reader.Read(&r)) != Reader::end; )
	{
		if (result == Reader::error)
		{
			cout << reader.message << endl;
			ok = false;
			break;
		}
		if (result == Reader::not_relation)
		{
			num_failed++;
			continue;
		}
		if (mode == "b" && (r.a < INT32_MIN || r.a > INT32_MAX))
		{
			num_too_large++;
			continue;
		}
		if (!write_relation(outfile, mode, side, r))
		{
			cout << "Could not write to the output file " << outputfile << "." << endl;
			ok = false;
			break;
		}
		num_relations++;
		if (num_relations % 10000000 == 0)
		{
			cout << "Converted " << num_relations << " relations." << endl;
		}
	}
	if (fclose(outfile) != 0 && ok)
	{
		cout << "Could not write to the output file " << outputfile << "." << endl;
		ok = false;
	}

	cout << "Converted " << num_relations << " relations (" << side_name(side) << " side)." << endl;
	if (reader.incomplete_last_line)
	{
		cout << "Skipped an incomplete last line of the input file." << endl;
	}
	if (num_failed > 0)
	{
		cout << "Skipped " << num_failed << " lines that are not relations (or have b of 2^32 or more, or more than " << MAX_PRIMES << " primes on the side)." << endl;
	}
	if (num_too_large > 0)
	{
		cout << "Skipped " << num_too_large << " relations with |a| of 2^31 or more, which the binary format cannot hold. Use the .ab format (mode a) for them." << endl;
	}
	return ok ? 0 : -1;
}
