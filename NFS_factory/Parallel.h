#pragma once
//A parallel loop whose threads are started for the loop and end with it. OpenMP's threads spin while they wait for the next
//loop, which in the batch pipeline, with short loops in several threads at once, took about as much CPU time as the loops.

#include <atomic>
#include <thread>
#include <vector>
#include <omp.h>

//body(i) for i from 0 to n - 1, by up to the number of threads of the threads option (omp_get_max_threads), which take chunk
//indices at a time, in order.
template <class Body> void parallel_for(int n, int chunk, Body body)
{
	const int chunks = (n + chunk - 1) / chunk;
	int threads = omp_get_max_threads();
	threads = (threads < chunks) ? threads : chunks;
	if (threads <= 1)
	{
		for (int i = 0; i < n; i++)
		{
			body(i);
		}
		return;
	}

	std::atomic<int> next(0);
	auto run = [&]()
	{
		for (;;)
		{
			const int start = next.fetch_add(chunk);
			if (start >= n)
			{
				return;
			}
			const int end = (n - start > chunk) ? start + chunk : n;
			for (int i = start; i < end; i++)
			{
				body(i);
			}
		}
	};
	std::vector<std::thread> pool;
	for (int k = 1; k < threads; k++)
	{
		pool.emplace_back(run);
	}
	run();
	for (std::thread& t : pool)
	{
		t.join();
	}
}
