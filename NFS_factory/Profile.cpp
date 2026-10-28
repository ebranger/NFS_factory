#include "Profile.h"
#include <iostream>
#include <iomanip>
#include <windows.h>
#include <psapi.h>

using namespace std;

unsigned long long profile_cycles[PROF_NUM_PHASES] = { 0 };

static const char* phase_names[PROF_NUM_PHASES] = {
	"Setup",
	"Reading relations",
	"Polynomial values",
	"Product tree",
	"Prime product mod root",
	"Scaling at root",
	"Remainder tree",
	"Leaves (square, gcd)",
	"Cofactorization",
	"Factoring relations",
	"Other side factors",
	"Writing relations"
};

//Print the time spent in each phase, as a share of the total, and the peak memory use.
void print_profile(unsigned long long total_cycles)
{
	unsigned long long accounted = 0;
	for (int i = 0; i < PROF_NUM_PHASES; i++)
	{
		accounted += profile_cycles[i];
	}

	cout << endl << "Time profile (Gcycles, share of total):" << endl;
	cout << fixed << setprecision(1);

	for (int i = 0; i < PROF_NUM_PHASES; i++)
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

	PROCESS_MEMORY_COUNTERS memory;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &memory, sizeof(memory)))
	{
		cout << "Peak memory use: " << memory.PeakWorkingSetSize / (1024 * 1024) << " MB" << endl;
	}

	cout.unsetf(ios::fixed);
	cout << setprecision(6);
}
