#include "Profile.h"
#include <iostream>
#include <iomanip>
#include <windows.h>
#include <psapi.h>

using namespace std;

profile_counter profile_cycles[PROF_NUM_PHASES] = {};
profile_counter profile_cycles_3lp[PROF_NUM_PHASES] = {};

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

//Print the time spent in each phase, as a share of the total, and the peak memory use. When the batches are checked in a
//pipeline, the phases overlap, and their shares can add up to more than 100 %.
void print_profile(unsigned long long total_cycles)
{
	unsigned long long accounted = 0;
	for (int i = 0; i < PROF_FIRST_WAIT; i++)
	{
		accounted += profile_cycles[i];
	}

	cout << endl << "Time profile (Gcycles, share of total):" << endl;
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
		cout << "3LP second stage (Gcycles, share of total):" << endl;
		for (int i = 0; i < PROF_FIRST_WAIT; i++)
		{
			if (profile_cycles_3lp[i] > 0)
			{
				cout << "  " << left << setw(25) << phase_names[i] << right << setw(8) << profile_cycles_3lp[i] / 1e9 << setw(7) << 100.0 * profile_cycles_3lp[i] / total_cycles << " %" << endl;
			}
		}
	}

	bool waits = false;
	for (int i = PROF_FIRST_WAIT; i < PROF_NUM_PHASES; i++)
	{
		waits = waits || profile_cycles[i] > 0;
	}
	if (waits)
	{
		cout << "Waits of the batch pipeline (Gcycles, share of total):" << endl;
		for (int i = PROF_FIRST_WAIT; i < PROF_NUM_PHASES; i++)
		{
			cout << "  " << left << setw(25) << phase_names[i] << right << setw(8) << profile_cycles[i] / 1e9 << setw(7) << 100.0 * profile_cycles[i] / total_cycles << " %" << endl;
		}
	}

	PROCESS_MEMORY_COUNTERS memory;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &memory, sizeof(memory)))
	{
		cout << "Peak memory use: " << memory.PeakWorkingSetSize / (1024 * 1024) << " MB" << endl;
	}

	cout.unsetf(ios::fixed);
	cout << setprecision(6);
}
