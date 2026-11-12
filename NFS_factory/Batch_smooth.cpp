#include "Batch_smooth.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <cerrno>
#include <time.h>
#include <limits.h>
#include "cofactorize.h"

#include "Profile.h"
#include "Prime_product.h"
#include "Smooth_factor.h"
#include "Parallel.h"

using namespace std;

#define MAX_TREE_LEVELS 20 //Default for the maximum number of levels, can be changed with treelevels in the polynomial file.
//max 22 requires ~3GB of RAM. However, for SNFS-220 I seem to get the same speed for lower size, which saves a lot of RAM. 20 is ~7% slower, but saved 2GB of RAM, so worth it.
//For larger factorizations, consider increasing this a little, but for every +1 the memory requirement for the tree roughly doubles.

#define MAX_CHARS_IN_RELATION_LIST 320 //Longest saved factor list. Relations with longer lists are rejected. Typical lists are less than 100 characters.

#define MAX_FACTORS_IN_RELATION 512 //Each factor is at least 2, so this is enough for any relation value of up to 512 bits.

static void append_prime(string& list, unsigned long long prime);

//Initialize the batch smoothness class
Batch_smooth::Batch_smooth(long long fblim, int lbp, int mfb, int max_bits_in_input, string filename, Polynomial* input_poly, int batch_checking_side, int print_limit, Polynomial* other_side_poly, bool rebuild_other_side, int max_tree_levels, tree_backend_kind backend, bool three_lp_stage, int pipeline)
{
	long long max_prime_limit = 1LL << lbp;

	this->three_lp_stage = three_lp_stage;
	if (three_lp_stage)
	{
		//The cofactors have no prime factors up to the fblim of the first stage, but the primes below that are cheap to include,
		//and then the candidates do not depend on the fblim used to find them.
		fblim = max_prime_limit - 1;
		mfb = 3 * lbp;
	}

	//if fblim and 2^lbp match, then we batch-check all numbers against the full primorial and do no cofactorization of any remaining numbers.
	// fblim > 2^lbp makes little sense? so assume equal, otherwise we will end up printing relations with some factors too large. Maybe check this in input file and warn user?
	if (max_prime_limit >= fblim && !three_lp_stage)
	{
		do_cofactorization = true;
	}
	else
	{ 
		do_cofactorization = false;
	}

	poly = input_poly;
	other_poly = other_side_poly;
	rebuild_all = rebuild_other_side && other_side_poly != NULL;
	profile = three_lp_stage ? profile_cycles_3lp : profile_cycles;

	if (poly->poly_is_binomial_SNFS_deg4())
	{
		using_snfs_deg4_binomial = true;
	}
	else
	{
		using_snfs_deg4_binomial = false;
	}

	//The coefficients of a binomial, for the native evaluation, if they fit in 63 bits.
	binomial_fits = using_snfs_deg4_binomial && mpz_cmpabs_ui(poly->get_coefficient(4).get_mpz_t(), (1ULL << 63) - 1) <= 0 && mpz_cmpabs_ui(poly->get_coefficient(0).get_mpz_t(), (1ULL << 63) - 1) <= 0;
	binomial_c4 = binomial_fits ? mpz_get_si(poly->get_coefficient(4).get_mpz_t()) : 0;
	binomial_c0 = binomial_fits ? mpz_get_si(poly->get_coefficient(0).get_mpz_t()) : 0;

	side = batch_checking_side;

	prime_product = new mpz_class();
	*prime_product = 1;

	Setup_prime_product(fblim, backend == TREE_BACKEND_GPU);
	max_batch_prime = fblim;

	//For handling 2lp cofactorization.
	one_lp_bit_limit = lbp;
	two_lp_min_limit = 2 * (long long)floor(log2((double)fblim));
	max_cofactor_bit_size = min(mfb, 2 * lbp);

	//Larger cofactors, up to mfb bits, can have three large primes. They must be composite and at least fblim^3.
	three_lp_max_bits = (mfb > 2 * lbp && !three_lp_stage) ? mfb : 0;
	mpz_ui_pow_ui(three_lp_min.get_mpz_t(), (mpir_ui)fblim, 3);
	second_stage = NULL;
	candidate_file = NULL;
	num_candidates = 0;
	num_candidates_skipped = 0;

	//Estiamte a good tree size for the remainder tree
	//Choose the smallest tree where the product of the relation values at the root is at least half the size of the prime product.
	//A larger tree makes the prime product modulo the root cheaper, but costs one more level of products and twice the memory.
	//For SNFS-240 on an i7-14700F, the root being 0.57 times the prime product (19 levels) was as fast as 1.15 times (20 levels),
	//and 0.29 times (18 levels) about 5 % slower.
	tree_size = (int)ceil(log2((double)mpz_sizeinbase(prime_product->get_mpz_t(), 2) / max_bits_in_input));
	const int tree_cap = (max_tree_levels > 0) ? max_tree_levels : MAX_TREE_LEVELS;
	//On the GPU, up to two levels more if that still leaves the most batches in progress (the memory of the GPU and of the computer
	//allows it): the root, whose Newton inverse is most of its cost, is then shared by more relations. For SNFS-240 on an RTX 4070
	//Ti SUPER, with alim 10^8, 20 levels (the root 1.15 times the prime product) were about 4 % faster than 19, with about 1.8 GB
	//per batch; with alim 2^24.7, the trees took 155, 129 and 130 ms per million relations with 18, 19 and 20 levels (the root 1.2,
	//2.4 and 4.8 times the prime product), and the smallest tree here has 17.
	for (int extra = 0; extra < 2; extra++)
	{
		if (backend == TREE_BACKEND_GPU && !three_lp_stage && pipeline <= 0 && tree_size + 1 <= tree_cap &&
			default_pipeline(backend, tree_size + 1, max_bits_in_input, false) >= GPU_MAX_PIPELINE)
		{
			tree_size++;
		}
	}
	tree_size = min(tree_cap, tree_size);  // Force the remainder tree to not be too big to save memory.
	tree_size = max(2, tree_size); // The remainder tree needs at least one level above the relations.

	max_entries_before_check = 1 << (tree_size - 1);
	save_filename = filename;
	output = NULL;
	num_relations_found = 0;
	num_relations_written = 0;
	num_relations_checked = 0;
	num_factoring_failures = 0;
	num_other_side_failures = 0;
	min_saved_prime_limit = print_limit; //Do not print primes less than this. Msieve reconstructs primes less than 1000, CADO needs all primes (so set this to 1 to print everything).

	factor1 = new mpz_class;
	factor2 = new mpz_class;

	//The batch slots, and the threads of the pipeline.
	num_relations_per_batch = max_entries_before_check;
	if (pipeline <= 0)
	{
		pipeline = three_lp_stage ? 1 : default_pipeline(backend, tree_size, max_bits_in_input);
	}
	pipelined = pipeline > 1 && !three_lp_stage;
	const int num_slots = pipelined ? pipeline : 1;
	for (int i = 0; i < num_slots; i++)
	{
		slots.push_back(New_slot(tree_size, max_bits_in_input, backend));
		free_slots.push_back(slots.back());
	}
	filling = free_slots.front();
	free_slots.pop_front();
	read_position[0] = read_position[1] = 0;
	submitted = 0;
	checked = 0;
	stopping = false;
	reader_blocked = 0;
	//The trees check the cofactors themselves, if they can (the GPU), with the limits of Classify_cofactor.
	if (do_cofactorization && mpz_size(three_lp_min.get_mpz_t()) <= 2)
	{
		Leaf_check check;
		check.max_cofactor_bits = max_cofactor_bit_size;
		check.three_lp_max_bits = three_lp_max_bits;
		check.two_lp_min_bits = (int)two_lp_min_limit;
		check.one_lp_bits = one_lp_bit_limit;
		check.three_lp_min[0] = mpz_getlimbn(three_lp_min.get_mpz_t(), 0);
		check.three_lp_min[1] = mpz_getlimbn(three_lp_min.get_mpz_t(), 1);
		check.max_batch_prime = (unsigned long long)max_batch_prime;
		for (Batch_slot* s : slots)
		{
			s->tree->Set_leaf_check(check);
		}
	}
	for (Batch_slot* s : slots)
	{
		s->tree->Set_leaf_factoring();
	}
	//The other side's norms on the GPU, if the backend can.
	other_factoring = NULL;
	if (other_poly != NULL)
	{
		vector<mpz_class> coefficients;
		for (int i = 0; i <= other_poly->get_degree(); i++)
		{
			coefficients.push_back(other_poly->get_coefficient((uint8)i));
		}
		other_factoring = create_other_side_factoring(backend, coefficients);
	}
	backend_kind = backend;
	writer = this;
	writer_busy = false;
	writer_stopping = false;
	flush_requested = false;
	output_size = -1;
	other_stage = NULL;
	other_stage_tried = other_stage_finished = false;
	other_tickets = 0;
	sink_next = 0;

	if (pipelined)
	{
		//The slots that are not being filled or checked can be in the trees at once (up to 4).
		const int tree_workers = min(4, max(1, num_slots - 2));
		for (int i = 0; i < tree_workers; i++)
		{
			workers.emplace_back(&Batch_smooth::Tree_worker, this);
		}
		workers.emplace_back(&Batch_smooth::Check_worker, this);
		writer_thread = thread(&Batch_smooth::Writer_worker, this);
	}
}

//Destructor, free all memory that was requested. Relations must already have been checked with Do_batch_check.
Batch_smooth::~Batch_smooth(void)
{
	{
		lock_guard<mutex> guard(lock);
		stopping = true;
	}
	changed.notify_all();
	for (thread& worker : workers)
	{
		worker.join();
	}
	//Only this object's own writer: that of a second stage is gone by now, and was flushed before.
	if (writer == this)
	{
		Flush_writer();
	}
	{
		lock_guard<mutex> guard(output_lock);
		writer_stopping = true;
	}
	output_changed.notify_all();
	if (writer_thread.joinable())
	{
		writer_thread.join();
	}

	if (num_factoring_failures > 0)
	{
		cout << "Could not fully factor " << num_factoring_failures << " smooth relations, they were not saved." << endl;
	}

	if (num_other_side_failures > 0)
	{
		cout << num_other_side_failures << " smooth relations had a factor list for the other side that did not match its norm, they were not saved." << endl;
	}

	delete second_stage;
	delete other_stage;
	delete other_factoring;
	if (candidate_file != NULL)
	{
		fclose(candidate_file);
	}
	if (output != NULL)
	{
		fclose(output);
	}

	for (Batch_slot* s : slots)
	{
		Delete_slot(s);
	}
	delete prime_product;
	delete factor1;
	delete factor2;
}

//Setup the primorials to be used.
void Batch_smooth::Setup_prime_product(long long smooth_bound, bool use_gpu)
{
	//So apparently mpir knows primorials. Problem solved!

	if (smooth_bound < 1000)
	{
		//I expect a minimum bound of ~2^20, anything less is probably a bad idea. 
		//Use 1000 as a basic sanity check. 
		return;
	}

	//Only the primes that can divide the norms are included, with small primes to higher powers.
	setup_prime_product(poly, (uint32)smooth_bound, prime_product, use_gpu);

	size_t size = mpz_sizeinbase(prime_product->get_mpz_t(), 2);
	console_message("Size of prime product: " + to_string(size) + " bits.");

	return;
}

//A batch slot, with the memory for the relations of a batch, and its tree, which allocates all its memory up front.
Batch_slot* Batch_smooth::New_slot(int tree_levels, int max_bits_in_input, tree_backend_kind backend)
{
	Batch_slot* s = new Batch_slot;
	s->tree = create_tree_backend(backend, tree_levels, max_bits_in_input);
	const int n = num_relations_per_batch;
	s->a_vals = new long long[n];
	s->b_vals = new long long[n];
	s->is_smooth = new bool[n];
	s->lp1 = new unsigned long long[n];
	s->lp2 = new unsigned long long[n];
	//The saved factor lists, in one block.
	s->list_of_factors = new char*[n];
	char* lists = new char[(size_t)n * MAX_CHARS_IN_RELATION_LIST];
	for (int i = 0; i < n; i++)
	{
		s->list_of_factors[i] = lists + (size_t)i * MAX_CHARS_IN_RELATION_LIST;
		s->a_vals[i] = 0;
		s->b_vals[i] = 0;
		s->lp1[i] = 0;
		s->lp2[i] = 0;
		s->is_smooth[i] = false;
	}
	s->binary_lists = false;
	s->entries = 0;
	s->sequence = 0;
	s->position[0] = s->position[1] = 0;
	s->root_bits = 0;
	return s;
}

void Batch_smooth::Delete_slot(Batch_slot* s)
{
	delete s->tree;
	delete[] s->list_of_factors[0];
	delete[] s->list_of_factors;
	delete[] s->a_vals;
	delete[] s->b_vals;
	delete[] s->is_smooth;
	delete[] s->lp1;
	delete[] s->lp2;
	delete s;
}

//Add one relation to be batch-smoothness checked. Returns false if the factor list is too long to be saved.
bool Batch_smooth::Add_number(long long a_val, long long b_val, const char* relation_factor_list)
{
	//add one relation for batch-factoring.
	if (strlen(relation_factor_list) >= MAX_CHARS_IN_RELATION_LIST)
	{
		return false;
	}

	Batch_slot& s = *filling;
	s.a_vals[s.entries] = a_val;
	s.b_vals[s.entries] = b_val;
	s.is_smooth[s.entries] = true; //Until the batch check shows otherwise.
	strcpy(s.list_of_factors[s.entries], relation_factor_list);
	s.binary_lists = false;

	s.entries++;

	//Is the tree full?
	if (s.entries == max_entries_before_check) //entries points at first empty position.
	{
		Submit();
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

	Batch_slot& s = *filling;
	s.a_vals[s.entries] = a_val;
	s.b_vals[s.entries] = b_val;
	s.is_smooth[s.entries] = true; //Until the batch check shows otherwise.

	//The saved list holds the number of primes, followed by the primes.
	unsigned int* saved = reinterpret_cast<unsigned int*>(s.list_of_factors[s.entries]);
	saved[0] = num_primes;
	if (num_primes > 0)
	{
		memcpy(saved + 1, primes, num_primes * sizeof(unsigned int));
	}
	s.binary_lists = true;

	s.entries++;

	//Is the tree full?
	if (s.entries == max_entries_before_check)
	{
		Submit();
	}
	return true;
}

//Add one 3LP candidate to the second stage: a relation with the factor list of the other side, and the cofactor of the batch checked side
//that is left after dividing out the primes up to fblim of the first stage. Returns false if there are too many primes to be saved.
bool Batch_smooth::Add_candidate(long long a_val, long long b_val, const unsigned int* primes, unsigned int num_primes, mpz_srcptr cofactor)
{
	if (num_primes > MAX_SAVED_PRIMES)
	{
		return false;
	}

	//The cofactor is checked instead of the polynomial value.
	mpz_set(filling->tree->Leaf(filling->entries)->get_mpz_t(), cofactor);
	return Add_number_binary(a_val, b_val, primes, num_primes);
}

bool Batch_smooth::Set_three_lp_handling(Batch_smooth* second_stage, string candidate_filename)
{
	this->second_stage = second_stage;
	if (second_stage != NULL)
	{
		second_stage->writer = this;
	}

	if (!candidate_filename.empty())
	{
		candidate_file = fopen(candidate_filename.c_str(), "ab");
		if (candidate_file == NULL)
		{
			return false;
		}
	}
	return true;
}

//Pass on the 3LP candidates of the batch that was just checked, to the second stage and/or the candidate file.
//The candidate file has the binary relation format with the cofactor added: (int)a, (uint)b, (uint)num_primes, (uint)prime1 ... (uint)prime_num_primes,
//(uint64)low and (uint64)high 64 bits of the cofactor.
void Batch_smooth::Send_candidates(Batch_slot& s)
{
	if (result_sink)
	{
		return; //an other_stage gives its candidates to the writer, in Output_solutions
	}
	vector<unsigned long long> primes;
	unsigned int saved[MAX_SAVED_PRIMES];
	long long* a_vals = s.a_vals;
	long long* b_vals = s.b_vals;

	for (size_t n = 0; n < s.candidate_index.size(); n++)
	{
		int index = s.candidate_index[n];
		Get_saved_primes(s, index, primes);

		bool fits = primes.size() <= MAX_SAVED_PRIMES && a_vals[index] >= INT_MIN && a_vals[index] <= INT_MAX && b_vals[index] <= UINT_MAX;
		for (size_t j = 0; j < primes.size() && fits; j++)
		{
			fits = primes[j] <= UINT_MAX;
			saved[j] = (unsigned int)primes[j];
		}

		if (!fits)
		{
			num_candidates_skipped++;
			continue;
		}

		num_candidates++;
		unsigned int num_primes = (unsigned int)primes.size();

		if (candidate_file != NULL)
		{
			int a = (int)a_vals[index];
			unsigned int b = (unsigned int)b_vals[index];
			unsigned long long cofactor[2] = { mpz_getlimbn(s.candidate_cofactor[n].get_mpz_t(), 0), mpz_getlimbn(s.candidate_cofactor[n].get_mpz_t(), 1) };
			fwrite(&a, sizeof(int), 1, candidate_file);
			fwrite(&b, sizeof(unsigned int), 1, candidate_file);
			fwrite(&num_primes, sizeof(unsigned int), 1, candidate_file);
			fwrite(saved, sizeof(unsigned int), num_primes, candidate_file);
			fwrite(cofactor, sizeof(unsigned long long), 2, candidate_file);
		}

		if (second_stage != NULL)
		{
			second_stage->Add_candidate(a_vals[index], b_vals[index], saved, num_primes, s.candidate_cofactor[n].get_mpz_t());
		}
	}

	s.candidate_index.clear();
}

void Batch_smooth::Finish()
{
	//The last batch, even if it is not full, and then all batches done.
	if (filling->entries > 0)
	{
		Submit();
	}
	if (pipelined)
	{
		unsigned long long t = profile_now();
		unique_lock<mutex> guard(lock);
		changed.wait(guard, [this] { return checked == submitted; });
		reader_blocked += profile_add(PROF_WAIT_READER, t) - t;
	}

	if (three_lp_max_bits > 0 && !result_sink)
	{
		string text = "Found " + to_string(num_candidates) + " 3LP candidates";
		if (num_candidates_skipped > 0)
		{
			text += " (and skipped " + to_string(num_candidates_skipped) + " with a, b or saved primes too large for the candidate format)";
		}
		console_message(text + ".");
	}

	if (second_stage != NULL)
	{
		second_stage->Finish();
	}
	Flush_writer();
	if (second_stage != NULL)
	{
		console_message("The 3LP second stage found " + to_string(second_stage->Get_num_relations_found()) + " relations.");
	}
}

//The saved factor list of the other side of a relation, as text.
void Batch_smooth::Get_saved_list(Batch_slot& s, int index, string& list)
{
	if (!s.binary_lists)
	{
		list = s.list_of_factors[index];
		return;
	}

	const unsigned int* saved = reinterpret_cast<const unsigned int*>(s.list_of_factors[index]);
	list.clear();
	for (unsigned int i = 1; i <= saved[0]; i++)
	{
		append_prime(list, saved[i]);
	}
}

//Whether a relation was added from a binary file without a factor list of the other side (a file of only a and b). An empty list
//in a text file is printed as it is: such files often hold the primes of one side only.
bool Batch_smooth::Saved_list_missing(Batch_slot& s, int index)
{
	return s.binary_lists && reinterpret_cast<const unsigned int*>(s.list_of_factors[index])[0] == 0;
}

//The saved factor list of the other side of a relation, as numbers.
void Batch_smooth::Get_saved_primes(Batch_slot& s, int index, vector<unsigned long long>& primes)
{
	primes.clear();

	if (s.binary_lists)
	{
		const unsigned int* saved = reinterpret_cast<const unsigned int*>(s.list_of_factors[index]);
		for (unsigned int i = 1; i <= saved[0]; i++)
		{
			primes.push_back(saved[i]);
		}
		return;
	}

	const char* text = s.list_of_factors[index];
	while (*text != '\0')
	{
		char* end;
		primes.push_back(strtoull(text, &end, 16));
		text = (*end == ',') ? end + 1 : end;
	}
}

//Check whether a remainder factors into two numbers less than the lp limit.
//If it does, it is a valid relation!
//The large primes of the cofactor of a leaf, which is what is left of the relation value after dividing out the primes of the batch
//check: 1 (smooth), one large prime, two large primes (split with Pollard rho), or a 3LP candidate. Sets lp1 and lp2 of the
//relation. Thread safe: cofactors that 64-bit rho cannot split are left for Split_with_tinyqs.
Batch_smooth::cofactor_result Batch_smooth::Classify_cofactor(Batch_slot& s, mpz_class& number, unsigned int index)
{
	int cofactor_bits = (int)mpz_sizeinbase(number.get_mpz_t(), 2);

	if (cofactor_bits > max_cofactor_bit_size)
	{
		//Too big to be split into two primes both less than lp_limit. It can still have three, if it is composite and at least fblim^3.
		if (cofactor_bits <= three_lp_max_bits && number >= three_lp_min && !is_probable_prime(number.get_mpz_t()))
		{
			return COFACTOR_CANDIDATE;
		}
		return COFACTOR_NOT_SMOOTH;
	}

	if (cofactor_bits > two_lp_min_limit)
	{
		//Can potentially split into two large primes, small enough to be useful
		if (is_probable_prime(number.get_mpz_t()))
		{
			return COFACTOR_NOT_SMOOTH;
		}

		if (cofactor_bits <= 63 && mpz_odd_p(number.get_mpz_t()))
		{
			//Pollard rho with 64-bit arithmetic.
			unsigned long long n = mpz_get_ui(number.get_mpz_t()), factor;
			if (split_composite64(n, &factor))
			{
				return Check_two_primes(s, factor, n / factor, index);
			}
		}

		//either rho failed or the cofactor is too big for 64-bit rho, so do qs.
		return COFACTOR_TINYQS;
	}

	if (cofactor_bits > one_lp_bit_limit)
	{
		//Should be a prime and too big.
		return COFACTOR_NOT_SMOOTH;
	}

	if (number == 1)
	{
		//no large primes found, finished!
		return COFACTOR_SMOOTH;
	}

	//should be one large prime by now. Make sure it is actually a large prime as well.
	if (number > (int64)max_batch_prime && is_probable_prime(number.get_mpz_t()))
	{
		s.lp1[index] = mpz_get_ui(number.get_mpz_t());
		s.lp2[index] = 0;
	}
	else
	{
		//I have seen this happen when the remainder was a product of one large prime and one very small that occured to a high power in the relation,
		//and the small prime was not found by the batch code due to the high power.
		//Since this seems to happend very rarely, for now just do a basic check if the sizes of the two factors makes the relation valid,
		//and let the relation factoring code handle it later, without any lp hints from here.
		unsigned int factors[MAX_FACTORS_IN_RELATION];
		int count = factor_general(factors, number.get_mpz_t(), MAX_FACTORS_IN_RELATION);

		if (count < 0)
		{
			return COFACTOR_NOT_SMOOTH;
		}

		for (int i = 0; i < count; i++)
		{
			if ((unsigned long long)factors[i] >> one_lp_bit_limit != 0)
			{
				//factor is 2^lbp or larger.
				return COFACTOR_NOT_SMOOTH;
			}
		}
	}

	//if the remainder after dividing out all small primes is less than the limit, it is prime and smooth (or possibly composite with high-power primes in relation, but still OK.)
	return COFACTOR_SMOOTH;
}

//The two factors of a 2LP cofactor must be large primes: below the large prime bound, and prime. If it turns out to be a small prime
//that occured to a high power, it can cause problems later when it is assumed to be a prime larger than the rlim/alim.
Batch_smooth::cofactor_result Batch_smooth::Check_two_primes(Batch_slot& s, unsigned long long p, unsigned long long q, unsigned int index)
{
	mpz_class f1, f2;
	mpz_set_ui(f1.get_mpz_t(), (mpir_ui)p);
	mpz_set_ui(f2.get_mpz_t(), (mpir_ui)q);
	if ((int)mpz_sizeinbase(f1.get_mpz_t(), 2) > one_lp_bit_limit || (int)mpz_sizeinbase(f2.get_mpz_t(), 2) > one_lp_bit_limit)
	{
		return COFACTOR_NOT_SMOOTH;
	}
	if (!is_probable_prime(f1.get_mpz_t()) || !is_probable_prime(f2.get_mpz_t()))
	{
		return COFACTOR_NOT_SMOOTH;
	}

	//save the two found factors, the smaller first.
	s.lp1[index] = min(p, q);
	s.lp2[index] = max(p, q);
	return COFACTOR_SMOOTH;
}

//Split a 2LP cofactor that 64-bit rho could not with tinyqs, and check the two factors. Not thread safe.
bool Batch_smooth::Split_with_tinyqs(Batch_slot& s, mpz_class& number, unsigned int index)
{
	*factor1 = 0;
	*factor2 = 0;
	int qs_sucess = tinyqs(number.get_mpz_t(), factor1->get_mpz_t(), factor2->get_mpz_t());

	if (qs_sucess == 0)
	{
		//unsuccesful
		return false;
	}

	//check if the factors are too big.
	if (mpz_sizeinbase(factor1->get_mpz_t(), 2) > one_lp_bit_limit || mpz_sizeinbase(factor2->get_mpz_t(), 2) > one_lp_bit_limit)
	{
		return false;
	}
	return Check_two_primes(s, factor1->get_ui(), factor2->get_ui(), index) == COFACTOR_SMOOTH;
}

//|c| x^4 for |x| < 2^32 and |c| < 2^63, which is less than 2^191, as three limbs.
static void binomial_term(unsigned long long term[3], long long c, long long x)
{
	unsigned long long abs_x = (unsigned long long)(x < 0 ? -x : x);
	unsigned long long abs_c = (unsigned long long)(c < 0 ? -c : c);
	unsigned long long square = abs_x * abs_x;
	unsigned long long power_high;
	unsigned long long power_low = _umul128(square, square, &power_high);
	unsigned long long high0, high1;
	term[0] = _umul128(power_low, abs_c, &high0);
	unsigned long long middle = _umul128(power_high, abs_c, &high1);
	unsigned char carry = _addcarry_u64(0, middle, high0, &term[1]);
	_addcarry_u64(carry, high1, 0, &term[2]);
}

//The absolute value of the binomial c4 a^4 + c0 b^4 with native arithmetic. Returns false, and leaves the result alone, if the
//coefficients or a and b are too large for it.
bool Batch_smooth::Binomial_value_native(mpz_ptr result, long long a, long long b)
{
	const long long limit = 1LL << 32;
	if (!binomial_fits || a <= -limit || a >= limit || b <= -limit || b >= limit)
	{
		return false;
	}

	unsigned long long term4[3], term0[3], value[3];
	binomial_term(term4, binomial_c4, a);
	binomial_term(term0, binomial_c0, b);

	if ((binomial_c4 < 0) == (binomial_c0 < 0))
	{
		unsigned char carry = _addcarry_u64(0, term4[0], term0[0], &value[0]);
		carry = _addcarry_u64(carry, term4[1], term0[1], &value[1]);
		_addcarry_u64(carry, term4[2], term0[2], &value[2]);
	}
	else
	{
		//The difference of the terms, the smaller from the larger.
		unsigned long long* larger = term4;
		unsigned long long* smaller = term0;
		if (term0[2] > term4[2] || (term0[2] == term4[2] && (term0[1] > term4[1] || (term0[1] == term4[1] && term0[0] > term4[0]))))
		{
			larger = term0;
			smaller = term4;
		}
		unsigned char borrow = _subborrow_u64(0, larger[0], smaller[0], &value[0]);
		borrow = _subborrow_u64(borrow, larger[1], smaller[1], &value[1]);
		_subborrow_u64(borrow, larger[2], smaller[2], &value[2]);
	}

	if (result->_mp_alloc < 3)
	{
		mpz_realloc2(result, 3 * GMP_LIMB_BITS);
	}
	int size = 3;
	while (size > 0 && value[size - 1] == 0)
	{
		size--;
	}
	for (int i = 0; i < size; i++)
	{
		result->_mp_d[i] = value[i];
	}
	result->_mp_size = size;
	return true;
}

//Hand on the slot that was filled: check it now (pipeline 1), or pass it to the tree stage and take a free slot for the next
//relations, waiting for one if all are in use. The second stage counts all of its time as one phase of the profile, with its
//phases kept separately.
void Batch_smooth::Submit()
{
	Batch_slot* s = filling;
	s->position[0] = read_position[0];
	s->position[1] = read_position[1];
	unsigned long long t = profile_now();
	if (!pipelined)
	{
		s->sequence = submitted++;
		Run_trees(*s);
		Run_checks(*s);
		checked++;
		const unsigned long long now = three_lp_stage ? profile_add(PROF_3LP, t) : profile_now();
		reader_blocked += now - t;
		return;
	}

	unique_lock<mutex> guard(lock);
	s->sequence = submitted++;
	tree_queue.push_back(s);
	changed.notify_all();
	changed.wait(guard, [this] { return !free_slots.empty(); });
	filling = free_slots.front();
	free_slots.pop_front();
	guard.unlock();
	reader_blocked += profile_add(PROF_WAIT_READER, t) - t;
}

//The polynomial values of a batch (the second stage checks the cofactors that were added), the product tree of the values,
//the scaled remainder tree of the prime product P, and the cofactors of the values: the parts without the primes of P.
void Batch_smooth::Run_trees(Batch_slot& s)
{
	unsigned long long t = profile_now();

	for (int index = 0; index < s.entries && !three_lp_stage; index++)
	{
		mpz_class* value = s.tree->Leaf(index);

		//use different poly evaluation if polynomial has form k*x^4+l which is faster.
		if (using_snfs_deg4_binomial)
		{
			if (!Binomial_value_native(value->get_mpz_t(), s.a_vals[index], s.b_vals[index]))
			{
				poly->evaluate_binomial_polynomial_homogeneous_SNFS_deg4(value, s.a_vals[index], s.b_vals[index], &s.temp_number);
			}
		}
		else
		{
			poly->evaluate_polynomial_homogenous(value, (int64)s.a_vals[index], (int64)s.b_vals[index]);
		}

		if (*value < 0)
		{
			*value *= -1;
		}
	}
	t = profile_add(profile, PROF_POLY_EVAL, t);

	//Leaves without a relation (in a batch that is not full) must be 1.
	for (int index = s.entries; index < num_relations_per_batch; index++)
	{
		*s.tree->Leaf(index) = 1;
	}

	s.root_bits = s.tree->Build_product_tree();
	t = profile_add(profile, PROF_PRODUCT_TREE, t);

	s.tree->Scaled_root(prime_product->get_mpz_t());
	t = profile_add(profile, PROF_ROOT, t);

	s.tree->Remainder_tree();
	t = profile_add(profile, PROF_REMAINDER_TREE, t);

	s.tree->Leaf_cofactors();
	profile_add(profile, PROF_LEAVES, t);
}

//The checks of the cofactors of a batch, and its output: the cofactor N / gcd(N, (P mod N)^2) of each relation value N is 1 if N
//is smooth (without cofactorization), or is checked for large primes. The cofactors are checked in parallel (with the threads
//option), and the steps that are not thread safe afterwards, in leaf order. The batches are done one at a time, in their order.
void Batch_smooth::Run_checks(Batch_slot& s)
{
	unsigned long long t = profile_now();
	if (verbose_output)
	{
		console_message("Size of top node of product tree: " + to_string(s.root_bits) + " bits.");
	}

	s.leaf_result.resize(num_relations_per_batch);
	const unsigned char* classes = s.tree->Leaf_classes();
	if (classes != NULL)
	{
		//Checked by the tree: the large primes of the smooth ones, and the rest that the CPU checks.
		const unsigned long long* large_primes = s.tree->Leaf_large_primes();
		parallel_for(num_relations_per_batch, 4096, [&](int index)
		{
			if (classes[index] == LEAF_CHECK_ON_CPU)
			{
				s.leaf_result[index] = (unsigned char)Classify_cofactor(s, *s.tree->Leaf_cofactor(index), index);
				return;
			}
			s.leaf_result[index] = classes[index];
			if (classes[index] == LEAF_SMOOTH)
			{
				s.lp1[index] = large_primes[2 * index];
				s.lp2[index] = large_primes[2 * index + 1];
			}
		});
	}
	else
	{
		parallel_for(num_relations_per_batch, 256, [&](int index)
		{
			mpz_class& Q = *s.tree->Leaf_cofactor(index);

			if (do_cofactorization == true)
			{
				s.leaf_result[index] = (unsigned char)Classify_cofactor(s, Q, index);
			}
			else
			{
				//Do not do 2lp testing. If the cofactor is 1 then the relation is smooth
				s.leaf_result[index] = (unsigned char)((Q > 1) ? COFACTOR_NOT_SMOOTH : COFACTOR_SMOOTH);
			}
		});
	}

	for (int index = 0; index < num_relations_per_batch; index++)
	{
		switch (s.leaf_result[index])
		{
		case COFACTOR_SMOOTH:
			break;
		case COFACTOR_CANDIDATE:
			if (s.candidate_cofactor.size() <= s.candidate_index.size())
			{
				s.candidate_cofactor.emplace_back();
			}
			s.candidate_cofactor[s.candidate_index.size()] = *s.tree->Leaf_cofactor(index);
			s.candidate_index.push_back(index);
			s.is_smooth[index] = false;
			break;
		case COFACTOR_TINYQS:
			if (!Split_with_tinyqs(s, *s.tree->Leaf_cofactor(index), index))
			{
				s.is_smooth[index] = false;
			}
			break;
		default:
			s.is_smooth[index] = false;
			break;
		}
	}
	profile_add(profile, PROF_COFACTOR, t);

	Send_candidates(s);
	Output_solutions(s);
}

//A thread of the tree stage: the batches in the order they were submitted.
void Batch_smooth::Tree_worker()
{
	unique_lock<mutex> guard(lock);
	for (;;)
	{
		unsigned long long t = profile_now();
		changed.wait(guard, [this] { return stopping || !tree_queue.empty(); });
		if (tree_queue.empty())
		{
			return;
		}
		Batch_slot* s = tree_queue.front();
		tree_queue.pop_front();
		guard.unlock();
		profile_add(PROF_WAIT_TREES, t);

		Run_trees(*s);

		guard.lock();
		trees_done.push_back(s);
		changed.notify_all();
	}
}

//The thread of the checks: the batches whose trees are done, in their order, after which their slot is free.
void Batch_smooth::Check_worker()
{
	unique_lock<mutex> guard(lock);
	for (;;)
	{
		unsigned long long t = profile_now();
		Batch_slot* s = NULL;
		changed.wait(guard, [&]
		{
			for (Batch_slot* x : trees_done)
			{
				if (x->sequence == checked)
				{
					s = x;
					return true;
				}
			}
			return stopping;
		});
		if (s == NULL)
		{
			return;
		}
		trees_done.erase(find(trees_done.begin(), trees_done.end(), s));
		guard.unlock();
		profile_add(PROF_WAIT_CHECKS, t);

		Run_checks(*s);

		guard.lock();
		checked++;
		free_slots.push_back(s);
		changed.notify_all();
	}
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

//Rebuild the factor list of the side that was not batch checked from its norm, with all primes of at least min_saved_prime_limit and
//their multiplicities. primes: the saved list, if any, which is checked against the norm; primes left out by the siever (typically
//small ones) are added. factors: all prime factors of the norm, in increasing order, if they were found (on the GPU), and otherwise
//(num_factors -1) the norm is factored here. Returns false if a saved prime does not divide the norm, or if the norm could not be
//fully factored. Thread safe.
bool Batch_smooth::Complete_other_side(long long a, long long b, vector<unsigned long long>& primes, const unsigned int* factors, int num_factors, string& list)
{
	if (num_factors >= 0)
	{
		//All primes of the norm: the saved ones must be among them, as often as they are listed.
		sort(primes.begin(), primes.end());
		vector<unsigned long long> all(factors, factors + num_factors);
		if (!includes(all.begin(), all.end(), primes.begin(), primes.end()))
		{
			return false;
		}
		primes.swap(all);
	}
	else
	{
		mpz_class norm;
		other_poly->evaluate_polynomial_homogenous(&norm, (int64)a, (int64)b);
		mpz_abs(norm.get_mpz_t(), norm.get_mpz_t());

		if (norm == 0)
		{
			return false;
		}

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
			unsigned int other_side_factor_list[MAX_FACTORS_IN_RELATION];
			int num_other_factors = factor_general(other_side_factor_list, norm.get_mpz_t(), MAX_FACTORS_IN_RELATION);

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
	}

	list.clear();
	for (size_t j = 0; j < primes.size(); j++)
	{
		if (primes[j] >= (unsigned long long)min_saved_prime_limit)
		{
			append_prime(list, primes[j]);
		}
	}
	return true;
}

//Whether the primes of the saved factor list of the other side divide its norm, as often as they are listed.
bool Batch_smooth::Saved_list_divides(Batch_slot& s, int index)
{
	mpz_class norm;
	other_poly->evaluate_polynomial_homogenous(&norm, (int64)s.a_vals[index], (int64)s.b_vals[index]);
	if (norm == 0)
	{
		return false;
	}
	vector<unsigned long long> primes;
	Get_saved_primes(s, index, primes);
	for (size_t j = 0; j < primes.size(); j++)
	{
		if (primes[j] < 2 || !mpz_divisible_ui_p(norm.get_mpz_t(), (mpir_ui)primes[j]))
		{
			return false;
		}
		mpz_divexact_ui(norm.get_mpz_t(), norm.get_mpz_t(), (mpir_ui)primes[j]);
	}
	return true;
}

//The output line of a relation found smooth: the checked side is factored (its large primes are known), and the factor list of
//the other side is the saved one, or rebuilt from its norm. Thread safe. If that norm is factored in the background
//(other_factoring), waiting is set, and relation_text is only the factor list of the checked side, for the writer.
Batch_smooth::relation_status Batch_smooth::Format_relation(Batch_slot& s, int index, string& relation_text, bool& waiting)
{
	const long long* a_vals = s.a_vals;
	const long long* b_vals = s.b_vals;
	unsigned int large[MAX_COFACTOR_PRIMES];
	unsigned int factors[MAX_FACTORS_IN_RELATION];
	int num_factors = 0;
	int num_large = 0;
	string checked_side_list; //factors of the relation value that was batch checked.
	string other_side_list; //factors of the other side, as saved or rebuilt.

	if (keys_only)
	{
		//Smoothness is settled by the batch check and the cofactor split; only the pair is written.
		waiting = false;
		relation_text = to_string(a_vals[index]) + "," + to_string(b_vals[index]) + "\n";
		return RELATION_WRITTEN;
	}

	if (Factor_checked_side(s, index, large, num_large, factors, num_factors) != RELATION_WRITTEN)
	{
		return RELATION_FACTORING_FAILED;
	}

	waiting = false;
	if (other_poly != NULL && (rebuild_all || Saved_list_missing(s, index)))
	{
		if (other_factoring != NULL)
		{
			waiting = true;
		}
		else
		{
			vector<unsigned long long> primes;
			Get_saved_primes(s, index, primes);
			if (!Complete_other_side(a_vals[index], b_vals[index], primes, NULL, -1, other_side_list))
			{
				return RELATION_OTHER_SIDE_FAILED;
			}
		}
	}
	else
	{
		//The saved list must divide the norm: a relation whose a or b was stored wrongly (such as a of 2^31 or more in the int32
		//of a binary file) was checked with the wrong value, and is left out.
		if (other_poly != NULL && !Saved_list_divides(s, index))
		{
			return RELATION_OTHER_SIDE_FAILED;
		}
		Get_saved_list(s, index, other_side_list);
	}

	//The factor list for the batch checked side: large primes first, then the other factors.
	for (int j = 0; j < num_large; j++)
	{
		if (large[j] >= (unsigned int)min_saved_prime_limit)
		{
			append_prime(checked_side_list, large[j]);
		}
	}

	for (int j = 0; j < num_factors; j++)
	{
		if (factors[j] >= (unsigned int)min_saved_prime_limit)
		{
			append_prime(checked_side_list, factors[j]);
		}
	}

	if (waiting)
	{
		relation_text.swap(checked_side_list);
	}
	else
	{
		relation_text = Relation_line(a_vals[index], b_vals[index], checked_side_list, other_side_list);
	}
	return RELATION_WRITTEN;
}

//The prime factors of the value of the batch checked side of a relation found smooth: num_large large primes, and the num_factors
//factors of the rest. Returns RELATION_FACTORING_FAILED if they could not be found. Thread safe.
Batch_smooth::relation_status Batch_smooth::Factor_checked_side(Batch_slot& s, int index, unsigned int* large, int& num_large, unsigned int* factors, int& num_factors)
{
	mpz_class value = *s.tree->Leaf(index);
	const long long* a_vals = s.a_vals;
	const long long* b_vals = s.b_vals;
	mpz_class temp;
	num_factors = 0;
	num_large = 0;

	//negative values are handled by msieve, ignore here.
	if (value < 0)
	{
		value = -value;
	}

	//factor the relation. First, divide out the large primes, which are printed first.
	if (three_lp_stage)
	{
		//The value is the cofactor from the first stage, and all its primes are below 2^lbp. Those are the large primes,
		//and the rest of the polynomial value is factored as for other relations. The tree may have factored the cofactor already
		//(as factor_smooth would, which takes values of up to two limbs).
		const unsigned int* tree_large = NULL;
		num_large = (s.tree_factors && mpz_size(value.get_mpz_t()) <= 2) ? s.tree->Leaf_factors(index, &tree_large) : -1;
		if (num_large >= 0 && num_large <= MAX_COFACTOR_PRIMES)
		{
			copy(tree_large, tree_large + num_large, large);
		}
		else
		{
			num_large = factor_smooth(large, value.get_mpz_t(), MAX_COFACTOR_PRIMES);
		}
		poly->evaluate_polynomial_homogenous(&temp, (int64)a_vals[index], (int64)b_vals[index]);
		mpz_abs(temp.get_mpz_t(), temp.get_mpz_t());

		if (num_large < 0 || temp == 0 || !mpz_divisible_p(temp.get_mpz_t(), value.get_mpz_t()))
		{
			num_factors = -1;
		}
		else
		{
			mpz_divexact(value.get_mpz_t(), temp.get_mpz_t(), value.get_mpz_t());
		}
	}
	else
	{
		if (s.lp1[index] > 1)
		{
			large[num_large++] = (unsigned int)s.lp1[index];
			mpz_divexact_ui(value.get_mpz_t(), value.get_mpz_t(), s.lp1[index]);
		}

		if (s.lp2[index] > 1)
		{
			large[num_large++] = (unsigned int)s.lp2[index];
			mpz_divexact_ui(value.get_mpz_t(), value.get_mpz_t(), s.lp2[index]);
		}
	}

	//The factors of the rest, which the tree may have found already (as factor_general would). Otherwise Pollard rho: since the
	//rest is smooth, it is a decent choice.
	const unsigned int* tree_factors = NULL;
	const int num_tree_factors = (s.tree_factors && !three_lp_stage && num_factors == 0) ? s.tree->Leaf_factors(index, &tree_factors) : -1;
	if (num_tree_factors >= 0)
	{
		num_factors = num_tree_factors;
		copy(tree_factors, tree_factors + num_factors, factors);
	}
	else if (num_factors == 0 && value > 1)
	{
		num_factors = factor_general(factors, value.get_mpz_t(), MAX_FACTORS_IN_RELATION);
	}

	//Typically: one factor is too big for an uint32. An incomplete relation is useless, so skip it.
	return (num_factors < 0) ? RELATION_FACTORING_FAILED : RELATION_WRITTEN;
}

//The line of a relation: a,b:rational factors:algebraic factors
string Batch_smooth::Relation_line(long long a, long long b, const string& checked_side, const string& other_side)
{
	char ab_text[50];
	sprintf(ab_text, "%lld,%lld:", a, b);
	string line = ab_text;
	//With side 0, the rational factors are those of the checked side, otherwise the algebraic ones.
	line += (side == 0) ? checked_side : other_side;
	line += ':';
	line += (side == 0) ? other_side : checked_side;
	line += '\n';
	return line;
}

//The lines of the relations of the batch that were found smooth, factored by the threads (with the threads option), in their order
//in the batch, for the writer. The norms of their other side that are factored in the background are submitted here. The slot is then
//cleared for the next batch.
void Batch_smooth::Output_solutions(Batch_slot& s)
{
	unsigned long long t = profile_now();

	if (result_sink)
	{
		//An other_stage: the factors of each relation, in the order they were added, for the writer that made it.
		s.tree_factors = s.tree->Leaf_factors_ready();
		vector<vector<unsigned long long>> found(s.entries);
		vector<unsigned char> factored(s.entries, 0);
		parallel_for(s.entries, 16, [&](int index)
		{
			unsigned int large[MAX_COFACTOR_PRIMES];
			unsigned int factors[MAX_FACTORS_IN_RELATION];
			int num_large = 0, num_factors = 0;
			if (s.is_smooth[index] && Factor_checked_side(s, index, large, num_large, factors, num_factors) == RELATION_WRITTEN)
			{
				found[index].assign(large, large + num_large);
				found[index].insert(found[index].end(), factors, factors + num_factors);
				sort(found[index].begin(), found[index].end());
				factored[index] = 1;
			}
		});
		vector<int> candidate(s.entries, -1);
		for (size_t n = 0; n < s.candidate_index.size(); n++)
		{
			candidate[s.candidate_index[n]] = (int)n;
		}
		for (int index = 0; index < s.entries; index++)
		{
			const bool is_candidate = candidate[index] >= 0;
			result_sink(sink_next + index, s.a_vals[index], s.b_vals[index], factored[index] != 0, found[index],
				is_candidate ? &s.candidate_cofactor[candidate[index]] : NULL, is_candidate ? s.tree->Leaf(index) : NULL);
			s.is_smooth[index] = false;
			s.lp1[index] = 0;
			s.lp2[index] = 0;
		}
		s.candidate_index.clear();
		sink_next += s.entries;
		s.entries = 0;
		profile_add(profile, PROF_FACTOR, t);
		return;
	}

	s.relation_texts.resize(num_relations_per_batch);
	s.relation_result.resize(num_relations_per_batch);
	s.relation_waiting.resize(num_relations_per_batch);
	s.tree_factors = s.tree->Leaf_factors_ready();
	parallel_for(s.entries, 16, [&](int index)
	{
		if (s.is_smooth[index])
		{
			bool waiting = false;
			s.relation_result[index] = (unsigned char)Format_relation(s, index, s.relation_texts[index], waiting);
			s.relation_waiting[index] = waiting;
		}
	});
	t = profile_add(profile, PROF_FACTOR, t);

	Output_batch* o = new Output_batch;
	o->source = this;
	for (int index = 0; index < s.entries; index++)
	{
		if (s.is_smooth[index])
		{
			if (s.relation_result[index] == RELATION_FACTORING_FAILED)
			{
				num_factoring_failures++;
			}
			else if (s.relation_result[index] == RELATION_OTHER_SIDE_FAILED)
			{
				num_other_side_failures++;
			}
			else if (s.relation_waiting[index])
			{
				Output_batch::Waiting w;
				w.line = o->lines.size();
				w.a = s.a_vals[index];
				w.b = s.b_vals[index];
				w.checked_side.swap(s.relation_texts[index]);
				Get_saved_primes(s, index, w.saved);
				w.ticket = -1;
				o->waiting.push_back(move(w));
				o->lines.emplace_back();
			}
			else
			{
				o->lines.push_back(move(s.relation_texts[index]));
			}
		}

		//reset for next iteration.
		s.is_smooth[index] = false;
		s.lp1[index] = 0;
		s.lp2[index] = 0;
	}
	o->batch_end = (writer == this);
	o->position[0] = s.position[0];
	o->position[1] = s.position[1];
	profile_add(profile, PROF_WRITE, t);

	if (verbose_output)
	{
		console_message(three_lp_stage ? "Ran a 3LP batch of " + to_string(s.entries) + " candidates, found " + to_string(o->lines.size()) + " relations" :
			"Ran a batch, found " + to_string(o->lines.size()) + " smooth numbers");
	}
	if (writer == this)
	{
		num_relations_checked += s.entries;
	}
	num_relations_found += (long long)o->lines.size();
	s.entries = 0;
	writer->Queue_output(o);
}

//Hand the relations found in a batch to the writer: its thread, or at once without a pipeline.
void Batch_smooth::Queue_output(Output_batch* o)
{
	if (!pipelined)
	{
		Writer_take(o);
		return;
	}
	lock_guard<mutex> guard(output_lock);
	output_queue.push_back(o);
	output_changed.notify_all();
}

//The thread of the writer: the batches of relations found, in order, and at the end (Flush_writer) the ones still held.
void Batch_smooth::Writer_worker()
{
	unique_lock<mutex> guard(output_lock);
	for (;;)
	{
		output_changed.wait(guard, [this] { return writer_stopping || flush_requested || !output_queue.empty(); });
		if (!output_queue.empty())
		{
			Output_batch* o = output_queue.front();
			output_queue.pop_front();
			writer_busy = true;
			guard.unlock();
			Writer_take(o);
			guard.lock();
			writer_busy = false;
			output_changed.notify_all();
			continue;
		}
		if (flush_requested)
		{
			writer_busy = true;
			guard.unlock();
			Write_held(true);
			guard.lock();
			writer_busy = false;
			flush_requested = false;
			output_changed.notify_all();
			continue;
		}
		return;
	}
}

//Write everything given to the writer, the batches still held included: the end of a run.
void Batch_smooth::Flush_writer()
{
	if (!writer->pipelined)
	{
		writer->Write_held(true);
		return;
	}
	unique_lock<mutex> guard(writer->output_lock);
	writer->flush_requested = true;
	writer->output_changed.notify_all();
	writer->output_changed.wait(guard, [this] { return writer->output_queue.empty() && !writer->writer_busy && !writer->flush_requested; });
}

//The other_stage of the writer, if the other side's norms are factored on the GPU.
void Batch_smooth::Create_other_stage()
{
	other_stage_tried = true;
	if (other_factoring == NULL || other_poly == NULL)
	{
		return;
	}
	//The bits of the norms: a and b of up to about 34 bits.
	int coefficient_bits = 0;
	for (int i = 0; i <= other_poly->get_degree(); i++)
	{
		coefficient_bits = max(coefficient_bits, (int)mpz_sizeinbase(other_poly->get_coefficient((uint8)i).get_mpz_t(), 2));
	}
	const int norm_bits = coefficient_bits + 34 * other_poly->get_degree() + 4;
	console_message("The factor lists of the other side of the relations found are batch checked against the primes up to 2^28.");
	other_stage = new Batch_smooth(OTHER_STAGE_BOUND, OTHER_STAGE_LARGE_PRIME_BITS, OTHER_STAGE_MAX_COFACTOR_BITS, norm_bits, save_filename, other_poly,
		1 - side, 1, NULL, false, OTHER_STAGE_LEVELS, backend_kind, false, 1);
	other_stage->profile = profile_cycles_other;
	other_stage->result_sink = [this](long long ticket, long long a, long long b, bool factored, vector<unsigned long long>& primes,
		const mpz_class* cofactor, const mpz_class* norm)
	{
		if (factored)
		{
			Other_result& r = other_results[ticket];
			r.state = OTHER_FACTORED;
			r.primes.swap(primes);
		}
		else if (cofactor != NULL && mpz_size(norm->get_mpz_t()) <= 4 && mpz_divisible_p(norm->get_mpz_t(), cofactor->get_mpz_t()))
		{
			//The rest of the norm and the cofactor, factored apart.
			mpz_class rest;
			mpz_divexact(rest.get_mpz_t(), norm->get_mpz_t(), cofactor->get_mpz_t());
			const mpz_class* parts[2] = { &rest, cofactor };
			for (const mpz_class* v : parts)
			{
				const size_t n = mpz_size(v->get_mpz_t());
				for (size_t k = 0; k < 4; k++)
				{
					split_limbs.push_back((k < n) ? mpz_getlimbn(v->get_mpz_t(), k) : 0);
				}
				split_sizes.push_back((unsigned int)n);
			}
			split_tickets.push_back(ticket);
		}
		else
		{
			fallback_tickets.push_back(ticket);
			fallback_a.push_back(a);
			fallback_b.push_back(b);
		}
	};
}

//The jobs of other_factoring for the norms that other_stage did not factor, after one of its trees: the cofactors and the rest of
//the norms apart, and the whole norms.
void Batch_smooth::Submit_fallback()
{
	if (!split_tickets.empty())
	{
		const int job = other_factoring->Submit_values(split_limbs.data(), split_sizes.data(), (int)split_sizes.size());
		for (long long ticket : split_tickets)
		{
			Other_result& r = other_results[ticket];
			r.state = OTHER_IN_JOB;
			r.job = job;
		}
		job_tickets[job].swap(split_tickets);
		job_split[job] = true;
		pending_jobs.push_back(job);
		split_tickets.clear();
		split_limbs.clear();
		split_sizes.clear();
	}
	if (fallback_tickets.empty())
	{
		return;
	}
	const int job = other_factoring->Submit(fallback_a.data(), fallback_b.data(), (int)fallback_a.size());
	for (long long ticket : fallback_tickets)
	{
		Other_result& r = other_results[ticket];
		r.state = OTHER_IN_JOB;
		r.job = job;
	}
	job_tickets[job].swap(fallback_tickets);
	job_split[job] = false;
	pending_jobs.push_back(job);
	fallback_tickets.clear();
	fallback_a.clear();
	fallback_b.clear();
}

//The results of a job of other_factoring, for all of its norms, after those of the jobs submitted before it: the results must be
//taken in that order, and a job of values and a job of norms from the same tree have tickets in between each other's.
void Batch_smooth::Take_job(int job)
{
	while (!pending_jobs.empty())
	{
		const int next = pending_jobs.front();
		pending_jobs.pop_front();
		Take_one_job(next);
		if (next == job)
		{
			return;
		}
	}
}

void Batch_smooth::Take_one_job(int job)
{
	vector<int> counts;
	vector<unsigned int> offsets, primes;
	other_factoring->Results(job, counts, offsets, primes);
	const vector<long long>& tickets = job_tickets[job];
	const bool split = job_split[job];
	for (size_t i = 0; i < tickets.size(); i++)
	{
		Other_result& r = other_results[tickets[i]];
		const size_t first = split ? 2 * i : i, values = split ? 2 : 1;
		r.state = OTHER_FACTORED;
		r.primes.clear();
		for (size_t v = first; v < first + values; v++)
		{
			if (counts[v] < 0)
			{
				r.state = OTHER_FAILED;
				break;
			}
			r.primes.insert(r.primes.end(), primes.begin() + offsets[v], primes.begin() + offsets[v] + counts[v]);
		}
		sort(r.primes.begin(), r.primes.end());
	}
	job_tickets.erase(job);
	job_split.erase(job);
}

//A batch of relations found, in the writer: the norms of the lines that wait for their other side go to other_stage, and the
//batches are written in order once the factors of all their lines are there.
void Batch_smooth::Writer_take(Output_batch* o)
{
	unsigned long long t = profile_now();
	if (!o->waiting.empty() && !other_stage_tried)
	{
		Create_other_stage();
	}
	for (Output_batch::Waiting& w : o->waiting)
	{
		if (other_stage != NULL)
		{
			w.ticket = other_tickets++;
			other_stage->Add_number_binary(w.a, w.b, NULL, 0);
			Submit_fallback();
		}
	}
	held.push_back(o);
	profile_add(profile, PROF_OTHER_SIDE, t);
	Write_held(false);
}

//Whether the factors of all the waiting lines of a batch are there.
bool Batch_smooth::Batch_ready(const Output_batch& o)
{
	for (const Output_batch::Waiting& w : o.waiting)
	{
		if (w.ticket >= 0 && other_results.find(w.ticket) == other_results.end())
		{
			return false;
		}
	}
	return true;
}

//Write the batches held, in order, as far as their factors are there; with flush, all of them, after the last norms are checked.
void Batch_smooth::Write_held(bool flush)
{
	if (flush && other_stage != NULL && !other_stage_finished)
	{
		unsigned long long t = profile_now();
		other_stage->Finish();
		other_stage_finished = true;
		Submit_fallback();
		profile_add(profile, PROF_OTHER_SIDE, t);
	}
	while (!held.empty() && (flush || Batch_ready(*held.front())))
	{
		Write_batch(*held.front());
		delete held.front();
		held.pop_front();
	}
}

//Finish the lines of a batch that waited for their other side, write them, and after a batch of this object, call the batch done
//callback with the size of the output file. The norms that other_stage did not factor are factored by other_factoring (on the GPU),
//or on the CPU if that fails too.
void Batch_smooth::Write_batch(Output_batch& o)
{
	unsigned long long t = profile_now();
	Batch_smooth& source = *o.source;
	vector<size_t> rest; //waiting lines without factors
	for (size_t k = 0; k < o.waiting.size(); k++)
	{
		Output_batch::Waiting& w = o.waiting[k];
		auto r = (w.ticket >= 0) ? other_results.find(w.ticket) : other_results.end();
		if (r != other_results.end() && r->second.state == OTHER_IN_JOB)
		{
			Take_job(r->second.job);
		}
		if (r == other_results.end())
		{
			rest.push_back(k);
		}
		else if (r->second.state == OTHER_FAILED)
		{
			//Not factored on the GPU at all: on the CPU (rare).
			string other_side;
			if (source.Complete_other_side(w.a, w.b, w.saved, NULL, -1, other_side))
			{
				o.lines[w.line] = source.Relation_line(w.a, w.b, w.checked_side, other_side);
			}
			else
			{
				source.num_other_side_failures++;
				source.num_relations_found--;
			}
		}
		else
		{
			vector<unsigned int> factors;
			for (unsigned long long p : r->second.primes)
			{
				factors.push_back((unsigned int)p);
			}
			string other_side;
			if (source.Complete_other_side(w.a, w.b, w.saved, factors.data(), (int)factors.size(), other_side))
			{
				o.lines[w.line] = source.Relation_line(w.a, w.b, w.checked_side, other_side);
			}
			else
			{
				source.num_other_side_failures++;
				source.num_relations_found--;
			}
		}
		if (r != other_results.end())
		{
			other_results.erase(r);
		}
	}
	if (!rest.empty())
	{
		vector<long long> a, b;
		for (size_t k : rest)
		{
			a.push_back(o.waiting[k].a);
			b.push_back(o.waiting[k].b);
		}
		vector<int> counts(rest.size(), -1);
		vector<unsigned int> offsets, primes;
		if (source.other_factoring != NULL)
		{
			source.other_factoring->Results(source.other_factoring->Submit(a.data(), b.data(), (int)a.size()), counts, offsets, primes);
		}
		for (size_t j = 0; j < rest.size(); j++)
		{
			Output_batch::Waiting& w = o.waiting[rest[j]];
			string other_side;
			if (source.Complete_other_side(w.a, w.b, w.saved, (counts[j] >= 0) ? primes.data() + offsets[j] : NULL, counts[j], other_side))
			{
				o.lines[w.line] = source.Relation_line(w.a, w.b, w.checked_side, other_side);
			}
			else
			{
				source.num_other_side_failures++;
				source.num_relations_found--;
			}
		}
	}
	t = profile_add(profile, PROF_OTHER_SIDE, t);

	if (output == NULL)
	{
		output = fopen(save_filename.c_str(), "a");
		if (output == NULL)
		{
			console_message("Could not open output file for writing relations.");
		}
	}

	int written = 0;
	output_size = -1;
	if (output != NULL)
	{
		for (const string& line : o.lines)
		{
			if (!line.empty())
			{
				fputs(line.c_str(), output);
				written++;
			}
		}
		fflush(output);
#ifdef _WIN32
		_fseeki64(output, 0, SEEK_END);
		output_size = _ftelli64(output);
#else
		fseeko(output, 0, SEEK_END);
		output_size = (long long)ftello(output);
#endif
	}
	num_relations_written += written;
	profile_add(profile, PROF_WRITE, t);

	if (o.batch_end && batch_done)
	{
		batch_done(o.position, output_size, num_relations_written);
	}
}

void Batch_smooth::Set_read_position(long long first, long long second)
{
	read_position[0] = first;
	read_position[1] = second;
}

void Batch_smooth::Set_batch_done(function<void(const long long position[2], long long output_size, long long written)> callback)
{
	batch_done = callback;
}

unsigned long long Batch_smooth::Get_reader_blocked_cycles()
{
	return reader_blocked;
}

long long Batch_smooth::Get_num_3lp_relations_found()
{
	return (second_stage != NULL) ? second_stage->Get_num_relations_found() : 0;
}

long long Batch_smooth::Get_num_relations_checked()
{
	return num_relations_checked;
}

long long Batch_smooth::Get_num_relations_found()
{
	return num_relations_found + ((second_stage != NULL) ? second_stage->Get_num_relations_found() : 0);
}