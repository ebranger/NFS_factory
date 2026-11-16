#include "Profile.h"
#include <iostream>
#include <iomanip>
#include <mutex>
#include <chrono>
#include <stdio.h>
#ifdef _WIN32
#include <io.h>
#include <windows.h>
#include <psapi.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#endif

using namespace std;

profile_counter profile_cycles[PROF_NUM_PHASES] = {};
profile_counter profile_cycles_3lp[PROF_NUM_PHASES] = {};
profile_counter profile_cycles_other[PROF_NUM_PHASES] = {};

static const char* phase_names[PROF_NUM_PHASES] = {
	"Setup",
	"Reading relations",
	"Polynomial values",
	"Product tree",
	"Root of remainder tree",
	"Remainder tree",
	"Leaves (square, gcd)",
	"Cofactorization",
	"Factoring relations",
	"Other side factors",
	"Writing relations",
	"3LP second stage",
	"Reader, for a slot",
	"Trees, for a batch",
	"Checks, for the trees"
};

//Print the time spent in each phase, as a share of the time of the whole run, and the peak memory use. The time of a phase is
//added up over the threads that do it: when the batches are checked in a pipeline, several threads work at once (up to four
//batches in the trees, the reader and the thread of the checks), so a phase can take more than 100 %, and the shares add up to
//more than 100 %.
void print_profile(unsigned long long total_cycles)
{
	unsigned long long accounted = 0;
	for (int i = 0; i < PROF_FIRST_WAIT; i++)
	{
		accounted += profile_cycles[i];
	}

	bool waits = false;
	for (int i = PROF_FIRST_WAIT; i < PROF_NUM_PHASES; i++)
	{
		waits = waits || profile_cycles[i] > 0;
	}

	cout << endl << "Time profile (Gcycles, share of the time of the run):" << endl;
	if (waits)
	{
		cout << "The time of each step is added up over the threads doing it at once, so a step can take more than 100 %." << endl;
	}
	cout << fixed << setprecision(1);

	for (int i = 0; i < PROF_FIRST_WAIT; i++)
	{
		if (profile_cycles[i] > 0)
		{
			cout << "  " << left << setw(25) << phase_names[i] << right << setw(8) << profile_cycles[i] / 1e9 << setw(7) << 100.0 * profile_cycles[i] / total_cycles << " %" << endl;
		}
	}

	if (total_cycles > accounted)
	{
		cout << "  " << left << setw(25) << "Other" << right << setw(8) << (total_cycles - accounted) / 1e9 << setw(7) << 100.0 * (total_cycles - accounted) / total_cycles << " %" << endl;
	}

	if (profile_cycles[PROF_3LP] > 0)
	{
		cout << "3LP second stage (Gcycles, share of the time of the run):" << endl;
		for (int i = 0; i < PROF_FIRST_WAIT; i++)
		{
			if (profile_cycles_3lp[i] > 0)
			{
				cout << "  " << left << setw(25) << phase_names[i] << right << setw(8) << profile_cycles_3lp[i] / 1e9 << setw(7) << 100.0 * profile_cycles_3lp[i] / total_cycles << " %" << endl;
			}
		}
	}

	bool other = false;
	for (int i = 0; i < PROF_FIRST_WAIT; i++)
	{
		other = other || profile_cycles_other[i] > 0;
	}
	if (other)
	{
		cout << "Batch check of the other side, in the writer (Gcycles, share of the time of the run):" << endl;
		for (int i = 0; i < PROF_FIRST_WAIT; i++)
		{
			if (profile_cycles_other[i] > 0)
			{
				cout << "  " << left << setw(25) << phase_names[i] << right << setw(8) << profile_cycles_other[i] / 1e9 << setw(7) << 100.0 * profile_cycles_other[i] / total_cycles << " %" << endl;
			}
		}
	}

	if (waits)
	{
		//The thread of the checks takes the batches in order, so it is the limit when it seldom waits for the trees. The trees (up
		//to four threads) are the limit when the reader waits for a slot and the checks wait for the trees.
		cout << "Waits of the batch pipeline (Gcycles, share of the time of the run):" << endl;
		for (int i = PROF_FIRST_WAIT; i < PROF_NUM_PHASES; i++)
		{
			cout << "  " << left << setw(25) << phase_names[i] << right << setw(8) << profile_cycles[i] / 1e9 << setw(7) << 100.0 * profile_cycles[i] / total_cycles << " %" << endl;
		}
	}

#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS memory;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &memory, sizeof(memory)))
	{
		cout << "Peak memory use: " << memory.PeakWorkingSetSize / (1024 * 1024) << " MB" << endl;
	}
#else
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage) == 0)
	{
		cout << "Peak memory use: " << usage.ru_maxrss / 1024 << " MB" << endl; //ru_maxrss is in kB on Linux
	}
#endif

	cout.unsetf(ios::fixed);
	cout << setprecision(6);
}

bool verbose_output = false;

static mutex console_lock;
static string status_text; //on the console
static size_t status_width = 0; //of the status line on the console, 0 if none is shown
static long long last_status = -1000000; //ms

static long long console_ms()
{
	return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now().time_since_epoch()).count();
}

bool console_is_terminal()
{
#ifdef _WIN32
	static const bool terminal = _isatty(_fileno(stdout)) != 0;
#else
	static const bool terminal = isatty(fileno(stdout)) != 0;
#endif
	return terminal;
}

//Rewrite the status line on the console. With console_lock.
static void draw_status()
{
	cout << '\r' << status_text;
	if (status_text.size() < status_width)
	{
		cout << string(status_width - status_text.size(), ' ');
	}
	cout << flush;
	status_width = (status_text.size() > status_width) ? status_text.size() : status_width;
}

void console_status(const string& text)
{
	lock_guard<mutex> guard(console_lock);
	const long long now = console_ms();
	if (now - last_status < (console_is_terminal() ? 500 : 60000))
	{
		return;
	}
	last_status = now;
	if (console_is_terminal())
	{
		status_text = text;
		draw_status();
	}
	else
	{
		cout << text << endl;
	}
}

void console_message(const string& text)
{
	lock_guard<mutex> guard(console_lock);
	if (status_width > 0)
	{
		//Over the status line, which is then written again below.
		cout << '\r' << text;
		if (text.size() < status_width)
		{
			cout << string(status_width - text.size(), ' ');
		}
		cout << endl;
		status_width = 0;
		draw_status();
	}
	else
	{
		cout << text << endl;
	}
}

void console_end_status(const string& text)
{
	lock_guard<mutex> guard(console_lock);
	if (console_is_terminal())
	{
		status_text = text;
		draw_status();
		cout << endl;
		status_width = 0;
		status_text.clear();
	}
	else
	{
		cout << text << endl;
	}
	last_status = console_ms();
}
