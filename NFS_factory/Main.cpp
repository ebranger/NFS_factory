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
#include <algorithm>
#include <thread>
#include <chrono>
#include <iomanip>
#include <omp.h>
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#define _fgetc_nolock getc_unlocked
#endif

using namespace std;

#define MAX_LINE_LENGTH 1024 //Longer lines in a relation file are skipped. gnfs-lasieve relations are typically less than 200 characters.
#define MAX_BINARY_PRIMES 32 //Maximum number of primes in one relation in a binary relation file. More than this is treated as a corrupt file.
#define READ_BUFFER_SIZE (16 * 1024 * 1024) //Binary relation files are read in blocks of this size.
#define READ_AHEAD_BLOCKS 4 //Reads of READ_BUFFER_SIZE in flight at once (Windows): a hard disk then reads at full speed.
#define AB_HEADER_BYTES 16 //The (a,b)-only relation format (.ab): see read_batch_from_ab_file.
#define AB_RECORD_BYTES 12
#define AB_VERSION 1
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
	bool keys_only; // "format: keys": write only a,b of the relations found (see Batch_smooth::keys_only).
	int tree_levels; // Maximum number of levels in the remainder tree. 0 means the default.
	int tree_levels_3lp; // Maximum number of levels in the remainder tree of a run on a 3LP candidate file.
	bool save_3lp; // 3LP candidates (mfb more than 2 lpb) are saved to a file instead of checked in a second stage during the run.
	string file_3lp; // File for the saved 3LP candidates. Empty means the default name.
	tree_backend_kind backend; // How the products of the trees are computed.
	int pipeline; // Batches in progress at once (see Batch_smooth). 0 means the default: on the GPU from its free memory, otherwise 1.
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
				else if (line == "keys")
				{
					params->cado_format = false;
					params->keys_only = true;
				}
				else
				{
					cout << "Unknown output format " << line << " on line " << line_number << " of the polynomial file. Use msieve, cado or keys." << endl;
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
			//Batches in progress at once: 1 checks each batch when it is full; more checks them in a pipeline of threads while the
			//next ones are read, each with its own tree.
			else if (data == "pipeline")
			{
				params->pipeline = stoi(line);
				if (params->pipeline < 1)
				{
					cout << "pipeline on line " << line_number << " of the polynomial file must be at least 1." << endl;
					return false;
				}
			}
			//How the products of the trees are computed: mpir (default of the CPU build), ntt, or gpu (default of the GPU build).
			else if (data == "backend")
			{
				if (line == "mpir")
				{
					params->backend = TREE_BACKEND_MPIR;
				}
				else if (line == "ntt")
				{
					params->backend = TREE_BACKEND_NTT;
				}
				else if (line == "gpu")
				{
#ifdef NFS_GPU
					params->backend = TREE_BACKEND_GPU;
#else
					cout << "This build has no GPU support, so backend gpu on line " << line_number << " of the polynomial file cannot be used. Use the GPU build (NFS_factory_gpu)." << endl;
					return false;
#endif
				}
				else
				{
					cout << "Unknown backend " << line << " on line " << line_number << " of the polynomial file. Use mpir, ntt or gpu." << endl;
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

//The size of a file in bytes, or -1 if it cannot be opened or is not a regular file (opening a named pipe here would take
//data from the program writing into it).
long long file_size(const string& name)
{
#ifdef _WIN32
	FILE* f = fopen(name.c_str(), "rb");
	if (f == NULL)
	{
		return -1;
	}
	_fseeki64(f, 0, SEEK_END);
	long long size = _ftelli64(f);
	fclose(f);
	return size;
#else
	struct stat st;
	if (stat(name.c_str(), &st) != 0 || !S_ISREG(st.st_mode))
	{
		return -1;
	}
	return (long long)st.st_size;
#endif
}

//Cut a file to size bytes. Returns false if that failed.
bool truncate_file(const string& name, long long size)
{
#ifdef _WIN32
	int fd = _open(name.c_str(), _O_RDWR | _O_BINARY);
	if (fd < 0)
	{
		return false;
	}
	bool ok = _chsize_s(fd, size) == 0;
	_close(fd);
	return ok;
#else
	return truncate(name.c_str(), (off_t)size) == 0;
#endif
}

//A number with its digits in groups of three, separated by commas.
string group_digits(long long value)
{
	string digits = to_string(value < 0 ? -value : value);
	string result;
	for (size_t i = 0; i < digits.size(); i++)
	{
		if (i > 0 && (digits.size() - i) % 3 == 0)
		{
			result += ',';
		}
		result += digits[i];
	}
	return (value < 0) ? "-" + result : result;
}

//The progress of a run over its relation files, for resuming an interrupted run and for stopping at a number of relations.
//After each batch, when its relations are checked and written to the output file, the position of its last relation is saved
//to the resume file: the input (the relation file, or the list of relation files), the file reached, and how many records (lines
//of a text file, relations of a binary file) of it are done, with the number of relations found so far. With a pipeline, the
//reader is ahead of that position, and the batches it read after it are checked again when the run resumes.
//The size of the output file is saved too: relations written after that point by an interrupted run are removed when it resumes,
//since their batch is checked again. 3LP candidates that wait for the second stage when a run is interrupted are not saved, and
//their relations are lost.
struct run_progress
{
	Batch_smooth* batch;
	string resume_filename;
	string output_filename;
	string input; //"file <name>" or "list <name>"
	vector<string> files;
	size_t file_index; //the file being read
	long long done_records; //records of that file that are done
	long long read_records; //records of that file read up to the last relation added
	long long found_before; //relations found before this run (when resuming)
	long long target; //stop when this many relations have been found (0: no target)
	long long relations_added; //by this run, for maxrels
	atomic<bool> target_reached; //set when the relations found reach the target (Check_target)
	bool complete; //all files were read

	long long Found()
	{
		return found_before + batch->Get_num_relations_found();
	}

	//The status line: the time since the reading started, the file of the last batch checked, the relations checked and found, and
	//their rates.
	long long start_ms;
	string Status(size_t file)
	{
		const double seconds = (double)(chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now().time_since_epoch()).count() - start_ms) / 1000;
		const long long elapsed = (long long)seconds;
		const long long checked = batch->Get_num_relations_checked();
		const long long found = batch->Get_num_relations_found();
		char text[512];
		string name = (file < files.size()) ? files[file] : files.back();
		const size_t slash = name.find_last_of("\\/");
		name = (slash == string::npos) ? name : name.substr(slash + 1);
		string where = name;
		if (files.size() > 1)
		{
			where += " (" + to_string(((file < files.size()) ? file : files.size() - 1) + 1) + " of " + to_string(files.size()) + ")";
		}
		const long long three_lp = batch->Get_num_3lp_relations_found();
		snprintf(text, sizeof(text), "%lld:%02lld:%02lld  %s  checked %.1fM, %.2fM/s  found %s%s, %.0f/s", elapsed / 3600, elapsed / 60 % 60, elapsed % 60,
			where.c_str(), checked / 1e6, (seconds > 0) ? checked / seconds / 1e6 : 0.0, group_digits(found_before + found).c_str(),
			(three_lp > 0) ? (" (3LP " + group_digits(three_lp) + ")").c_str() : "", (seconds > 0) ? found / seconds : 0.0);
		return text;
	}

	//The resume file is written by a thread of its own, since creating and renaming a file can take several ms on Windows. The
	//position and the size of the output file are taken together when a batch is done, and the latest of them is written; if the
	//run stops before that, the resume file holds an earlier one, which is just as valid.
	thread writer;
	mutex writer_lock;
	condition_variable writer_changed;
	string pending; //the next contents of the resume file, empty if none
	bool writing;
	bool writer_stop;

	void Start_writer()
	{
		writing = false;
		writer_stop = false;
		writer = thread([this]
		{
			unique_lock<mutex> guard(writer_lock);
			for (;;)
			{
				writer_changed.wait(guard, [this] { return writer_stop || !pending.empty(); });
				if (pending.empty())
				{
					return;
				}
				string text;
				text.swap(pending);
				writing = true;
				guard.unlock();
				Write(text);
				guard.lock();
				writing = false;
				writer_changed.notify_all();
			}
		});
	}

	//Wait until the resume file is written, and end the writer.
	void Stop_writer()
	{
		{
			lock_guard<mutex> guard(writer_lock);
			writer_stop = true;
		}
		writer_changed.notify_all();
		if (writer.joinable())
		{
			writer.join();
		}
	}

	void Save()
	{
		Save(file_index, done_records, file_size(output_filename), false, Found());
	}

	//The resume file, now or by the writer thread (background), with the size of the output file.
	void Save(size_t file, long long records, long long bytes, bool background, long long found)
	{
		char text[4096];
		snprintf(text, sizeof(text), "input %s\nfile %zu\nname %s\nrecords %lld\nfound %lld\noutput %lld\ncomplete %d\n", input.c_str(), file,
			(file < files.size()) ? files[file].c_str() : "-", records, found, (bytes > 0) ? bytes : 0, complete ? 1 : 0);
		if (background && writer.joinable())
		{
			{
				lock_guard<mutex> guard(writer_lock);
				pending = text;
			}
			writer_changed.notify_all();
			return;
		}
		//Written now: after any that the writer has.
		{
			unique_lock<mutex> guard(writer_lock);
			writer_changed.wait(guard, [this] { return !writing; });
			pending.clear();
		}
		Write(text);
	}

	void Write(const string& text)
	{
		string temp = resume_filename + ".tmp";
		FILE* f = fopen(temp.c_str(), "w");
		if (f == NULL)
		{
			cout << "Could not write the resume file " << temp << endl;
			return;
		}
		fputs(text.c_str(), f);
		fclose(f);
		remove(resume_filename.c_str());
		rename(temp.c_str(), resume_filename.c_str());
	}

	//Read a resume file of the same input. Returns false if there is none.
	bool Load()
	{
		ifstream f(resume_filename);
		string key, value, saved_input;
		size_t index = 0;
		long long records = 0, found = 0, bytes = 0;
		int done = 0;
		int keys = 0;
		while (f >> key && getline(f, value))
		{
			value = trim(value);
			try
			{
				if (key == "input") { saved_input = value; keys++; }
				else if (key == "file") { index = stoull(value); keys++; }
				else if (key == "records") { records = stoll(value); keys++; }
				else if (key == "found") { found = stoll(value); keys++; }
				else if (key == "output") { bytes = stoll(value); keys++; }
				else if (key == "complete") { done = stoi(value); keys++; }
			}
			catch (const exception&)
			{
				return false;
			}
		}
		if (keys != 6 || saved_input != input)
		{
			return false;
		}
		file_index = index;
		done_records = records;
		found_before = found;
		complete = (done != 0) || index >= files.size();

		//Relations written after the saved position are written again.
		long long size = file_size(output_filename);
		if (size > bytes && !truncate_file(output_filename, bytes))
		{
			cout << "Could not remove the relations written after the saved position from " << output_filename << endl;
		}
		if (size >= 0 && size < bytes)
		{
			cout << "The output file " << output_filename << " is shorter than when the resume file was written." << endl;
		}
		return true;
	}

	//Called after each batch, when its relations are written, with the position of its last relation (the file and the records
	//of it read) and the size of the output file.
	void Batch_done(const long long position[2], long long output_size, long long written)
	{
		Save((size_t)position[0], position[1], (output_size >= 0) ? output_size : file_size(output_filename), true, found_before + written);
		console_status(Status((size_t)position[0]));
		Check_target();
	}

	//Called before each relation is added to the batch, with the records of the current file read so far, for the position of
	//the relation.
	void Adding(long long records_read)
	{
		batch->Set_read_position((long long)file_index, records_read);
	}

	//Whether the relations found reach the target, which they count from when their batch is checked: when the other side is batch
	//checked, the relations are written up to a tree of that check later.
	bool Check_target()
	{
		if (target > 0 && Found() >= target)
		{
			if (!target_reached.exchange(true))
			{
				console_message("Found " + to_string(Found()) + " relations, at least the " + to_string(target) + " that were asked for.");
			}
			return true;
		}
		return false;
	}

	//Called after each relation added to the batch, with the records of the current file read so far. Returns true if
	//reading should stop.
	bool Relation_added(long long records_read, nfs_parameters* param)
	{
		relations_added++;
		read_records = records_read;
		if (target_reached || Check_target())
		{
			return true;
		}
		if (param->maxrels > 0 && relations_added >= param->maxrels)
		{
			console_message("Stopping after " + to_string(relations_added) + " relations (maxrels).");
			return true;
		}
		return false;
	}
};

enum read_result { READ_END, READ_STOPPED, READ_FAILED };

//Read the relations of a text relation file, from record skip_records on.
read_result read_batch_from_file(string filename, Batch_smooth* batch, nfs_parameters* param, run_progress& progress, long long skip_records)
{
	FILE* infile;

	infile = fopen(filename.c_str(), "r");

	if (infile == NULL)
	{
		console_message("Could not open input file " + filename);
		return READ_FAILED;
	}
	setvbuf(infile, NULL, _IOFBF, 1 << 20);

	char buffer[MAX_LINE_LENGTH];
	long long a, b;
	const char *tmp;
	char* ratlist, *alglist;

	long long relation_num = 0;
	long long failed_relations = 0;
	long long skipped_relations = 0;
	long long records = 0; //lines read
	bool incomplete_last_line = false;
	bool stopped = false;

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
			if (++records > skip_records)
			{
				failed_relations++;
			}
			continue;
		}

		//Lines done by an earlier run.
		if (++records <= skip_records)
		{
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
		progress.Adding(records);
		if (!batch->Add_number(a, b, (batch->side == 0) ? alglist : ratlist))
		{
			failed_relations++;
			continue;
		}

		relation_num++;

		if (progress.Relation_added(records, param))
		{
			stopped = true;
			break;
		}
	}

	fclose(infile);

	if (!stopped)
	{
		progress.done_records = records;
	}

	if (incomplete_last_line)
	{
		console_message("Skipped an incomplete last line in the relation file.");
	}

	if (failed_relations > 0)
	{
		console_message("Failed to parse " + to_string(failed_relations) + " relations out of a total of " + to_string(relation_num + failed_relations + skipped_relations) + " relations.");
	}

	if (skipped_relations > 0)
	{
		console_message("Skipped " + to_string(skipped_relations) + " relations with a or b larger than maxa or maxb.");
	}

	return stopped ? READ_STOPPED : READ_END;
}

//A file read in order, in blocks of READ_BUFFER_SIZE. On Windows, READ_AHEAD_BLOCKS reads are in flight at once (overlapped I/O),
//so that a hard disk keeps reading while the relations are parsed, or the reader waits for a batch: one read at a time, as fread
//does, got about 180 MB/s from a disk that reads 260 MB/s. Elsewhere, fread. A stream (Open_stream: standard input, or a pipe)
//is read with fread everywhere.
class Block_reader
{
public:
	Block_reader()
	{
		stream = NULL;
#ifdef _WIN32
		handle = INVALID_HANDLE_VALUE;
#else
		file = NULL;
#endif
	}
	~Block_reader()
	{
		Close();
	}

	//Read the file from byte offset on.
	bool Open(const string& name, long long offset = 0)
	{
#ifdef _WIN32
		handle = CreateFileA(name.c_str(), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_OVERLAPPED | FILE_FLAG_SEQUENTIAL_SCAN, NULL);
		if (handle == INVALID_HANDLE_VALUE)
		{
			return false;
		}
		next_offset = (unsigned long long)offset;
		current = 0;
		position = size = 0;
		finished = false;
		for (int i = 0; i < READ_AHEAD_BLOCKS; i++)
		{
			blocks[i].resize(READ_BUFFER_SIZE);
			ZeroMemory(&requests[i], sizeof(OVERLAPPED));
			requests[i].hEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
			Issue(i);
		}
		return true;
#else
		file = fopen(name.c_str(), "rb");
		if (file != NULL && offset > 0 && fseeko(file, (off_t)offset, SEEK_SET) != 0)
		{
			//Not seekable (a pipe): reading from the start instead would give the wrong relations.
			fclose(file);
			file = NULL;
		}
		return file != NULL;
#endif
	}

	//Read standard input ("-") or a file that is not a regular file (a named pipe) from its start, in order, with fread. It can
	//be opened only once: what was read cannot be read again, so the caller skips what it does not need by reading past it.
	bool Open_stream(const string& name)
	{
		if (name == "-")
		{
#ifdef _WIN32
			_setmode(_fileno(stdin), _O_BINARY);
#endif
			stream = stdin;
		}
		else
		{
			stream = fopen(name.c_str(), "rb");
		}
		return stream != NULL;
	}

	//Up to max of the next bytes of the file into out. Returns how many, 0 at the end of the file (or on a read error).
	size_t Read(unsigned char* out, size_t max)
	{
		if (stream != NULL)
		{
			return fread(out, 1, max, stream);
		}
#ifdef _WIN32
		size_t written = 0;
		while (written < max && !finished)
		{
			if (position == size)
			{
				//The next block, in the order the reads were issued.
				DWORD got = 0;
				if (!pending[current] || !GetOverlappedResult(handle, &requests[current], &got, TRUE))
				{
					got = 0;
				}
				pending[current] = false;
				position = 0;
				size = got;
				if (got == 0)
				{
					finished = true;
					break;
				}
			}
			const size_t n = std::min(max - written, size - position);
			memcpy(out + written, blocks[current].data() + position, n);
			written += n;
			position += n;
			if (position == size)
			{
				//The block is used up: read the next part of the file into it, unless this was the end.
				if (size == READ_BUFFER_SIZE)
				{
					Issue(current);
				}
				current = (current + 1) % READ_AHEAD_BLOCKS;
				if (size < READ_BUFFER_SIZE)
				{
					finished = true;
				}
				position = size = 0;
			}
		}
		return written;
#else
		return fread(out, 1, max, file);
#endif
	}

	void Close()
	{
		if (stream != NULL)
		{
			if (stream != stdin)
			{
				fclose(stream);
			}
			stream = NULL;
			return;
		}
#ifdef _WIN32
		if (handle == INVALID_HANDLE_VALUE)
		{
			return;
		}
		CancelIo(handle);
		for (int i = 0; i < READ_AHEAD_BLOCKS; i++)
		{
			DWORD got;
			if (pending[i])
			{
				GetOverlappedResult(handle, &requests[i], &got, TRUE);
			}
			CloseHandle(requests[i].hEvent);
		}
		CloseHandle(handle);
		handle = INVALID_HANDLE_VALUE;
#else
		if (file != NULL)
		{
			fclose(file);
			file = NULL;
		}
#endif
	}

private:
#ifdef _WIN32
	void Issue(int i)
	{
		requests[i].Offset = (DWORD)next_offset;
		requests[i].OffsetHigh = (DWORD)(next_offset >> 32);
		ResetEvent(requests[i].hEvent);
		pending[i] = ReadFile(handle, blocks[i].data(), READ_BUFFER_SIZE, NULL, &requests[i]) || GetLastError() == ERROR_IO_PENDING;
		next_offset += READ_BUFFER_SIZE;
	}

	HANDLE handle;
	vector<unsigned char> blocks[READ_AHEAD_BLOCKS];
	OVERLAPPED requests[READ_AHEAD_BLOCKS];
	bool pending[READ_AHEAD_BLOCKS];
	unsigned long long next_offset; //of the next read issued
	int current; //the block being used, in the order of the file
	size_t position, size; //in the current block
	bool finished;
#else
	FILE* file;
#endif
	FILE* stream; //standard input or a pipe (Open_stream), else NULL
};

//Whether a relation file is read as a stream (Block_reader::Open_stream): "-" is standard input, and a name that is not a regular
//file (a named pipe, or /dev/stdin) cannot be opened twice or seeked in either.
static bool is_stream_input(const string& name)
{
	if (name == "-")
	{
		return true;
	}
#ifdef _WIN32
	return false;
#else
	struct stat st;
	return stat(name.c_str(), &st) == 0 && !S_ISREG(st.st_mode);
#endif
}

//A buffer for reading a binary file in large blocks.
struct read_buffer
{
	Block_reader file;
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

			size_t read = file.Read(data.data() + filled, data.size() - filled);
			if (read == 0)
			{
				return false;
			}
			filled += read;
		}
		return true;
	}
};

//Read relations from a binary format file, from record skip_records on. With candidates, the file holds 3LP candidates, for a batch
//that is the 3LP second stage.
read_result read_batch_from_binary_file(string filename, Batch_smooth* batch, nfs_parameters* param, bool candidates, run_progress& progress, long long skip_records)
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
	long long records = 0; //relations read
	bool corrupt = false;
	bool stopped = false;

	read_buffer input;
	if (!input.file.Open(filename))
	{
		console_message("Failed to open input saved relations file " + filename);
		return READ_FAILED;
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

		//Relations done by an earlier run.
		if (++records <= skip_records)
		{
			continue;
		}

		if (!relation_in_range(a, b, param))
		{
			skipped_relations++;
			continue;
		}

		//Add the relation to be batch smoothness checked.
		progress.Adding(records);
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

		if (progress.Relation_added(records, param))
		{
			stopped = true;
			break;
		}
	}

	if (!corrupt && !stopped && input.filled > input.position)
	{
		//Part of a relation is left at the end of the file.
		corrupt = true;
	}
	input.file.Close();

	if (!stopped)
	{
		progress.done_records = records;
	}

	if (corrupt)
	{
		console_message("The binary relation file is truncated or corrupt after " + to_string(records) + " relations. The rest of the file was skipped.");
	}

	if (skipped_relations > 0)
	{
		console_message("Skipped " + to_string(skipped_relations) + " relations with a or b larger than maxa or maxb.");
	}

	return stopped ? READ_STOPPED : READ_END;
}

//Read relations from a file of only a and b (.ab), from record skip_records on. The file starts with a header of AB_HEADER_BYTES:
//"NFAB", (uint)version (AB_VERSION), (uint)the side whose primes the relations were sieved and saved for (0 rational, 1 algebraic,
//2 not known), and (uint)flags (bit 0: the relations were validated when the file was made). Then (int64)a, (uint32)b per relation,
//AB_RECORD_BYTES: a takes 64 bits, since the sievers' a can reach 2^31 and more, which the binary format cannot hold.
//The factor list of the side that is not batch checked is rebuilt from its norm for the relations found smooth (on the GPU in the
//GPU build), so the side in the header must not be the checked side.
read_result read_batch_from_ab_file(string filename, Batch_smooth* batch, nfs_parameters* param, run_progress& progress, long long skip_records)
{
	unsigned int header[AB_HEADER_BYTES / sizeof(unsigned int)];
	const bool stream = is_stream_input(filename);
	const string name = (filename == "-") ? "standard input" : filename;
	read_buffer input;
	input.data.resize(READ_BUFFER_SIZE);
	input.position = 0;
	input.filled = 0;
	bool header_read;
	if (stream)
	{
		//A stream (standard input, a pipe) is read once, from its start: the header, then the relations.
		if (!input.file.Open_stream(filename))
		{
			console_message("Failed to open input relation file " + name);
			return READ_FAILED;
		}
		header_read = input.ensure(AB_HEADER_BYTES);
		if (header_read)
		{
			memcpy(header, input.data.data(), AB_HEADER_BYTES);
			input.position = AB_HEADER_BYTES;
		}
	}
	else
	{
		FILE* f = fopen(filename.c_str(), "rb");
		if (f == NULL)
		{
			console_message("Failed to open input relation file " + filename);
			return READ_FAILED;
		}
		header_read = fread(header, AB_HEADER_BYTES, 1, f) == 1;
		fclose(f);
	}
	if (!header_read || memcmp(header, "NFAB", 4) != 0 || header[1] != AB_VERSION)
	{
		console_message(name + " is not a relation file of only a and b (.ab) of version " + to_string(AB_VERSION) + ".");
		return READ_FAILED;
	}
	if (header[2] == (unsigned int)param->smoothness_checking_side)
	{
		console_message("The relations of " + name + " were sieved for side " + to_string(header[2]) + ", which is the side to check here (side in the polynomial file).");
		return READ_FAILED;
	}

	long long relation_num = 0;
	long long skipped_relations = 0;
	long long records = skip_records; //relations read, those done by an earlier run included
	bool stopped = false;

	if (stream)
	{
		//The relations done by an earlier run are skipped by reading past them.
		unsigned long long skip = (unsigned long long)AB_RECORD_BYTES * (unsigned long long)skip_records;
		while (skip > 0 && input.ensure(1))
		{
			const size_t n = (size_t)std::min<unsigned long long>(skip, input.filled - input.position);
			input.position += n;
			skip -= n;
		}
		if (skip > 0)
		{
			console_message(name + " ended before the " + to_string(skip_records) + " relations done by the earlier run: it is not the same input.");
			return READ_FAILED;
		}
	}
	else if (!input.file.Open(filename, AB_HEADER_BYTES + AB_RECORD_BYTES * skip_records))
	{
		//The relations done by an earlier run are skipped by starting after them.
		console_message("Failed to open input relation file " + filename);
		return READ_FAILED;
	}

	while (input.ensure(AB_RECORD_BYTES))
	{
		long long a;
		unsigned int b;
		const unsigned char* record = input.data.data() + input.position;
		memcpy(&a, record, sizeof(long long));
		memcpy(&b, record + 8, sizeof(unsigned int));
		input.position += AB_RECORD_BYTES;
		records++;

		if (!relation_in_range(a, b, param))
		{
			skipped_relations++;
			continue;
		}

		progress.Adding(records);
		batch->Add_number_binary(a, b, NULL, 0);
		relation_num++;

		if (progress.Relation_added(records, param))
		{
			stopped = true;
			break;
		}
	}

	const bool partial = !stopped && input.filled > input.position;
	input.file.Close();

	if (!stopped)
	{
		progress.done_records = records;
	}

	if (partial)
	{
		console_message("The relation file " + name + " ends with part of a relation, which was skipped.");
	}

	if (skipped_relations > 0)
	{
		console_message("Skipped " + to_string(skipped_relations) + " relations with a or b larger than maxa or maxb.");
	}

	return stopped ? READ_STOPPED : READ_END;
}

//Whether name ends with suffix.
static bool ends_with(const string& name, const char* suffix)
{
	const size_t n = strlen(suffix);
	return name.size() >= n && name.compare(name.size() - n, n, suffix) == 0;
}

int main(int argc, char* argv[])
{
	//NFS_factory <polynomial file> <relation file> <output file> [-rels N]
	//NFS_factory <polynomial file> -list <list file> <output file> [-rels N]
	string inputfile, outputfile, polyfile, listfile;
	long long target_relations = 0;
#ifdef NFS_GPU
	//The GPU build runs as one process: 8 threads keep up with the GPU on SNFS-240 (RTX 4070 Ti SUPER), or fewer on a smaller CPU.
	int threads = (int)std::min(8u, std::max(1u, thread::hardware_concurrency()));
#else
	int threads = 1; //The CPU build is typically run as several processes.
#endif
	int gpu = 0;
	vector<string> positional;
	bool arguments_ok = true;
	for (int i = 1; i < argc; i++)
	{
		string arg = argv[i];
		if (arg == "-list" && i + 1 < argc)
		{
			listfile = argv[++i];
		}
		else if (arg == "-rels" && i + 1 < argc)
		{
			try
			{
				target_relations = stoll(argv[++i]);
			}
			catch (const exception&)
			{
				arguments_ok = false;
			}
		}
		else if (arg == "-gpu" && i + 1 < argc)
		{
			try
			{
				gpu = stoi(argv[++i]);
			}
			catch (const exception&)
			{
				arguments_ok = false;
			}
		}
		else if (arg == "-verbose")
		{
			verbose_output = true;
		}
		else if (arg == "-threads" && i + 1 < argc)
		{
			try
			{
				threads = stoi(argv[++i]);
			}
			catch (const exception&)
			{
				arguments_ok = false;
			}
			arguments_ok = arguments_ok && threads >= 1;
		}
		else
		{
			positional.push_back(arg);
		}
	}
	if (!arguments_ok || positional.size() != (listfile.empty() ? 3u : 2u))
	{
		cout << "Usage: NFS_factory <polynomial file> <relation file> <output file> [-rels N] [-threads N] [-gpu N] [-verbose]" << endl;
		cout << "   or: NFS_factory <polynomial file> -list <list file> <output file> [-rels N] [-threads N] [-gpu N] [-verbose]" << endl;
		cout << "The list file names relation files, one per line, which are read in order. With -rels, the run stops when N relations" << endl;
		cout << "have been found. An interrupted run continues where it was, from <output file>.resume, when it is run again with the" << endl;
		cout << "same relation file or list. A summary is written to <output file>.status at the end. -threads sets the number of" << endl;
		cout << "threads for checking the cofactors and factoring the relations found (default 1 in the CPU build, 8 in the GPU build)." << endl;
		cout << "-gpu selects the GPU of the GPU build, numbered from 0 (default 0). During the run, a status line shows the progress;" << endl;
		cout << "-verbose also prints a line for each batch." << endl;
		cout << "Relation files: text (a,b:rational primes:algebraic primes), binary (.bin), or only a and b (.ab), whose factor list of the" << endl;
		cout << "side that is not checked is rebuilt for the relations found (on the GPU in the GPU build). A .3lp file holds 3LP candidates." << endl;
		cout << "A relation file named - is an .ab file read from standard input, as in: xz -dc store.ab.xz | NFS_factory poly - output" << endl;
		return 0;
	}
	omp_set_num_threads(threads);

	polyfile = positional[0];
	inputfile = listfile.empty() ? positional[1] : listfile;
	outputfile = positional.back();

	run_progress progress;
	progress.input = (listfile.empty() ? "file " : "list ") + inputfile;
	if (listfile.empty())
	{
		progress.files.push_back(inputfile);
	}
	else
	{
		ifstream list(listfile);
		if (!list)
		{
			cout << "Could not open the list file " << listfile << endl;
			return -1;
		}
		string line;
		while (getline(list, line))
		{
			line = trim(line);
			if (!line.empty() && line[0] != '#')
			{
				progress.files.push_back(line);
			}
		}
		if (progress.files.empty())
		{
			cout << "The list file " << listfile << " names no relation files." << endl;
			return -1;
		}
	}
	if (count(progress.files.begin(), progress.files.end(), string("-")) > 1)
	{
		cout << "Standard input (-) can be read only once, but the list names it more than once." << endl;
		return -1;
	}
	progress.resume_filename = outputfile + ".resume";
	progress.output_filename = outputfile;
	progress.file_index = 0;
	progress.done_records = 0;
	progress.read_records = 0;
	progress.found_before = 0;
	progress.target = target_relations;
	progress.relations_added = 0;
	progress.target_reached = false;
	progress.complete = false;

	Batch_smooth* batch;

	//Wall-clock time: clock() is wall time with MSVC, but the CPU time of all threads on Linux.
	std::chrono::steady_clock::time_point time;
	unsigned long long program_start = profile_now();

	nfs_parameters* param = new nfs_parameters();
	param->maxa = LLONG_MAX; //No limit unless set in the polynomial file.
	param->maxb = LLONG_MAX;
	param->maxrels = 0;

	param->print_limit = 1000; //Msieve default.
	param->cado_format = false;
	param->keys_only = false;
	param->tree_levels = 0;
	param->tree_levels_3lp = STAGE2_ONLY_TREE_LEVELS;
	param->save_3lp = false;
	param->pipeline = 0;
#ifdef NFS_GPU
	param->backend = TREE_BACKEND_GPU;
#else
	param->backend = TREE_BACKEND_MPIR;
#endif

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

#ifdef NFS_GPU
	//The GPU, before anything else uses it.
	if (param->backend == TREE_BACKEND_GPU && !gpu_init(gpu))
	{
		return -1;
	}
#endif

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
	bool candidate_input = listfile.empty() && inputfile.size() > 4 && inputfile.compare(inputfile.size() - 4, 4, ".3lp") == 0;
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
		//CADO-NFS needs all primes on both sides, so print everything and rebuild the factor list of the other side of every relation.
		//Otherwise it is rebuilt only for the relations read without one (.ab files).
		param->print_limit = 1;
		cout << "Writing relations in CADO-NFS format, with all primes on both sides." << endl;
	}

	time = std::chrono::steady_clock::now();
	unsigned long long t = profile_now();
	if (candidate_input)
	{
		//The cofactors are less than 2^(3 lbp). The tree size is limited by 3lptreelevels instead of treelevels.
		cout << "Checking 3LP candidates against all primes up to 2^" << lbp << "." << endl;
		batch = new Batch_smooth(lim, lbp, mfb, 3 * lbp, outputfile, checked_poly, param->smoothness_checking_side, param->print_limit, other_poly, param->cado_format, param->tree_levels_3lp, param->backend, true);
		batch->keys_only = param->keys_only;
	}
	else
	{
		batch = new Batch_smooth(lim, lbp, mfb, tree_bits, outputfile, checked_poly, param->smoothness_checking_side, param->print_limit, other_poly, param->cado_format, param->tree_levels, param->backend, false, param->pipeline);
		batch->keys_only = param->keys_only;
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
			second_stage = new Batch_smooth(lim, lbp, mfb, 3 * lbp, outputfile, checked_poly, param->smoothness_checking_side, param->print_limit, other_poly, param->cado_format, param->tree_levels, param->backend, true);
			second_stage->keys_only = param->keys_only;
		}

		if (!batch->Set_three_lp_handling(second_stage, candidate_filename))
		{
			cout << "Could not open the 3LP candidate file " << candidate_filename << endl;
			delete batch;
			return -1;
		}
	}
	profile_add(PROF_SETUP, t);
	cout << "Prime products set up, it took " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - time).count() << " seconds." << endl;

	//What the run does: the batches, and the input.
	long long input_bytes = 0;
	for (const string& file : progress.files)
	{
		input_bytes += max(0LL, file_size(file));
	}
	cout << "Batches of " << group_digits(batch->Get_batch_size()) << " relations (" << batch->Get_tree_levels() << " tree levels), " << batch->Get_pipeline()
		<< ((batch->Get_pipeline() > 1) ? " in progress at once, " : " at a time, ") << threads << ((threads > 1) ? " threads." : " thread.") << endl;
	cout << "Input: " << progress.files.size() << ((progress.files.size() > 1) ? " relation files, " : " relation file, ") << fixed << setprecision(2)
		<< input_bytes / 1e9 << " GB. Output: " << outputfile << endl;
	cout.unsetf(ios::fixed);
	cout << setprecision(6);

	time = std::chrono::steady_clock::now();

	long long total_num_relations = 0;

	//The time in the readers that is not spent on batch checking, or waiting for it, is counted as reading.
	t = profile_now();

	progress.batch = batch;
	batch->Set_batch_done([&progress](const long long position[2], long long output_size, long long written) { progress.Batch_done(position, output_size, written); });
	progress.Start_writer();
	if (progress.Load())
	{
		if (progress.complete)
		{
			cout << "The resume file " << progress.resume_filename << " shows that all relation files were done, with " << progress.found_before << " relations found. Remove it to run again." << endl;
		}
		else
		{
			cout << "Resuming from " << progress.resume_filename << ": file " << progress.file_index + 1 << " of " << progress.files.size() << " (" << progress.files[progress.file_index]
				<< ") after " << progress.done_records << " records, with " << progress.found_before << " relations found before." << endl;
		}
	}
	else
	{
		progress.Save();
	}

	progress.start_ms = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now().time_since_epoch()).count();
	bool read_failed = false;
	bool stopped = progress.complete || (progress.target > 0 && progress.found_before >= progress.target);
	long long skip_records = progress.done_records;
	for (; !stopped && progress.file_index < progress.files.size(); progress.file_index++)
	{
		const string& file = progress.files[progress.file_index];
		read_result result;
		if (progress.files.size() > 1 && (verbose_output || !console_is_terminal()))
		{
			//On a console, the status line shows the file.
			console_message("Reading relation file " + to_string(progress.file_index + 1) + " of " + to_string(progress.files.size()) + ": " + file);
		}
		if (candidate_input)
		{
			//3LP candidates saved by an earlier run, with a .3lp extension.
			result = read_batch_from_binary_file(file, batch, param, true, progress, skip_records);
		}
		else if (ends_with(file, ".ab") || file == "-")
		{
			//Only a and b, with a .ab extension, or from standard input.
			result = read_batch_from_ab_file(file, batch, param, progress, skip_records);
		}
		else if (file.find(".bin") == string::npos)
		{
			//for reading gnfs-lasieve output.
			result = read_batch_from_file(file, batch, param, progress, skip_records);
		}
		else
		{
			//for reading binary formated files, assumed to have a .bin extension.
			result = read_batch_from_binary_file(file, batch, param, false, progress, skip_records);
		}
		skip_records = 0;

		if (result == READ_FAILED)
		{
			read_failed = true;
			break;
		}
		if (result == READ_STOPPED)
		{
			//After Finish below, the relations read so far are done.
			progress.done_records = progress.read_records;
			stopped = true;
			break;
		}
		//The relations of the next file continue the same batch.
		progress.done_records = 0;
	}

	//Check the remaining relations, even if the tree is not filled, and the 3LP candidates. Then all relations read are done.
	if (!read_failed)
	{
		batch->Finish();
		progress.complete = progress.complete || !stopped;
		if (progress.complete)
		{
			progress.file_index = progress.files.size();
			progress.done_records = 0;
		}
		progress.Save();
	}
	total_num_relations = batch->Get_num_relations_found();
	console_end_status(progress.Status(min(progress.file_index, progress.files.size() - 1)));

	unsigned long long reader_cycles = profile_now() - t;
	profile_cycles[PROF_READ] += reader_cycles - batch->Get_reader_blocked_cycles();

	cout << "Batch-checking all numbers took " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - time).count() << " seconds. Found a total of " << total_num_relations << " relations." << endl;

	//A summary for a driver script: relations found by this run and in total, whether the target was reached, and whether all files are done.
	string status_filename = outputfile + ".status";
	FILE* status = fopen(status_filename.c_str(), "w");
	if (status != NULL)
	{
		fprintf(status, "found %lld\ntotal %lld\ntarget %lld\ntarget_reached %d\nfiles_done %zu\nfiles %zu\ncomplete %d\nerror %d\n", total_num_relations,
			progress.Found(), progress.target, (progress.target > 0 && progress.Found() >= progress.target) ? 1 : 0,
			progress.complete ? progress.files.size() : progress.file_index, progress.files.size(), progress.complete ? 1 : 0, read_failed ? 1 : 0);
		fclose(status);
	}

	progress.Stop_writer();
	delete batch;
	delete param->algebraic_polynomial;
	delete param->rational_polynomial;
	delete param;

	if (read_failed)
	{
		//A relation file could not be read. The resume file shows where the run stopped.
		return -1;
	}

	print_profile(profile_now() - program_start);
	cout << "Batch checking complete." << endl;
	return 0;
}
