//The integer NTT on the GPU, see ntt_gpu.h.
//The transform kernels are those of gpuecm2's ntt.cu (by the same author). Lengths up to 4096 are done by one kernel (kRow), with
//each thread holding 8 elements and doing up to 3 radix-2 levels in registers per pass, exchanging data through shared memory
//between passes. Longer transforms use the four-step method T = R C: R-point transforms down the columns (kCol), a twist by
//omega_T^(l k1), then C-point transforms along the rows (kRow).

#include "ntt_gpu.h"
#include "arith32.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>

void cuda_check(cudaError_t error, const char* call, const char* file, int line)
{
	if (error != cudaSuccess)
	{
		fprintf(stderr, "CUDA error %s in %s at %s:%d\n", cudaGetErrorString(error), call, file, line);
		exit(3);
	}
}

namespace
{

constexpr int NTT_THREADS = 512;
constexpr int TILE = 4096; //elements per block

template <bool INV>
__device__ __forceinline__ void bfly(uint32_t& u, uint32_t& v, uint32_t w, uint32_t wp, uint32_t p)
{
	if (!INV)
	{
		uint32_t s = add_mod(u, v, p), d = sub_mod(u, v, p);
		u = s;
		v = mul_shoup(d, w, wp, p);
	}
	else
	{
		uint32_t t = mul_shoup(v, w, wp, p);
		v = sub_mod(u, t, p);
		u = add_mod(u, t, p);
	}
}

//A pass does R radix-2 levels (strides sb 2^(R-1) .. sb) of an M-point transform. Groups of G = 2^R elements {g + sb i}; each
//thread holds 8 / G groups. Element (u, i) of thread qq is at index elem_index(qq, u, i).
template <int R>
__device__ __forceinline__ int elem_index(int qq, int u, int i, int sb)
{
	constexpr int G = 1 << R, NG = 8 / G;
	int t = qq * NG + u;
	return (t / sb) * G * sb + (t % sb) + sb * i;
}

template <int R, bool INV>
__device__ __forceinline__ void do_pass(uint32_t v[8], int qq, int sb, int M, const uint32_t* sw, const uint32_t* swp, uint32_t p)
{
	constexpr int G = 1 << R, NG = 8 / G;
#pragma unroll
	for (int u = 0; u < NG; u++)
	{
		const int rr = (qq * NG + u) % sb;
		uint32_t* e = v + u * G;
#pragma unroll
		for (int step = 0; step < R; step++)
		{
			const int l = INV ? step : R - 1 - step; //level: half size s = sb << l
			const int f = M / (2 * (sb << l));
#pragma unroll
			for (int i = 0; i < G; i++)
			{
				if (i & (1 << l))
				{
					continue;
				}
				const int ex = (rr + sb * (i & ((1 << l) - 1))) * f;
				bfly<INV>(e[i], e[i + (1 << l)], sw[ex], swp[ex], p);
			}
		}
	}
}

template <int R>
__device__ __forceinline__ void index_set(int idx[8], int qq, int sb)
{
#pragma unroll
	for (int u = 0; u < 8 / (1 << R); u++)
	{
#pragma unroll
		for (int i = 0; i < (1 << R); i++)
		{
			idx[u * (1 << R) + i] = elem_index<R>(qq, u, i, sb);
		}
	}
}

//The pass schedule: forward = top levels first; the first forward pass (last inverse pass) has r0 = logM - 3 (npass - 1) levels.
__device__ __forceinline__ void pass_params(int logM, int k, bool inv, int& r, int& sb)
{
	const int npass = (logM + 2) / 3, r0 = logM - 3 * (npass - 1);
	const int kk = inv ? npass - 1 - k : k; //position in forward order
	r = (kk == 0) ? r0 : 3;
	const int levels_above = (kk == 0) ? 0 : r0 + 3 * (kk - 1); //levels done before this pass (forward)
	sb = 1 << (logM - levels_above - r);
}

template <bool INV>
__device__ __forceinline__ void run_pass(uint32_t v[8], int r, int qq, int sb, int M, const uint32_t* sw, const uint32_t* swp, uint32_t p)
{
	if (r == 3)
	{
		do_pass<3, INV>(v, qq, sb, M, sw, swp, p);
	}
	else if (r == 2)
	{
		do_pass<2, INV>(v, qq, sb, M, sw, swp, p);
	}
	else
	{
		do_pass<1, INV>(v, qq, sb, M, sw, swp, p);
	}
}

__device__ __forceinline__ void get_index_set(int idx[8], int r, int qq, int sb)
{
	if (r == 3)
	{
		index_set<3>(idx, qq, sb);
	}
	else if (r == 2)
	{
		index_set<2>(idx, qq, sb);
	}
	else
	{
		index_set<1>(idx, qq, sb);
	}
}

//M-point transforms of contiguous rows: X[(q nrows + row) M + i].
template <bool INV>
__global__ void __launch_bounds__(NTT_THREADS) kRow(uint32_t* __restrict__ X, size_t nrows, int logM, Gpu_ntt_tables t)
{
	extern __shared__ uint32_t sm[];
	const int M = 1 << logM, tpr = M / 8, rows_per_block = TILE / M;
	uint32_t* sw = sm;
	uint32_t* swp = sm + M / 2;
	uint32_t* sh = sm + M;
	const int q = blockIdx.y;
	const uint32_t p = t.p[q];
	const uint32_t* wt = INV ? t.iw : t.w;
	const uint32_t* wtp = INV ? t.iwp : t.wp;
	const int stride = 4096 / M;
	for (int i = threadIdx.x; i < M / 2; i += NTT_THREADS)
	{
		sw[i] = wt[(size_t)q * 2048 + i * stride];
		swp[i] = wtp[(size_t)q * 2048 + i * stride];
	}
	const int lrow = threadIdx.x / tpr, qq = threadIdx.x % tpr;
	const size_t row = (size_t)blockIdx.x * rows_per_block + lrow;
	const bool active = row < nrows;
	uint32_t* g = X + ((size_t)q * nrows + row) * M;
	uint32_t* s = sh + lrow * (M + M / 8);
	const int npass = (logM + 2) / 3;
	uint32_t v[8];
	int idx[8], r, sb;
	pass_params(logM, 0, INV, r, sb);
	get_index_set(idx, r, qq, sb);
	if (active)
	{
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			v[i] = g[idx[i]];
		}
	}
	__syncthreads();
	for (int k = 0; k < npass; k++)
	{
		if (active)
		{
			run_pass<INV>(v, r, qq, sb, M, sw, swp, p);
		}
		if (k == npass - 1)
		{
			break;
		}
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			s[idx[i] + (idx[i] >> 3)] = v[i];
		}
		__syncthreads();
		pass_params(logM, k + 1, INV, r, sb);
		get_index_set(idx, r, qq, sb);
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			v[i] = s[idx[i] + (idx[i] >> 3)];
		}
		__syncthreads();
	}
	if (active)
	{
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			g[idx[i]] = v[i];
		}
	}
}

//R-point transforms down the columns of each T = R C polynomial, with the four-step twist (after the forward pass, before the
//inverse pass). COLS = 4096 / R columns per block, R / 8 threads per column.
template <bool INV>
__global__ void __launch_bounds__(NTT_THREADS) kCol(uint32_t* __restrict__ X, size_t rows, int logT, int logR, Gpu_ntt_tables t)
{
	extern __shared__ uint32_t sm[];
	const int R = 1 << logR, COLS = TILE / R, logC = logT - logR;
	const size_t C = (size_t)1 << logC, T = (size_t)1 << logT;
	uint32_t* sw = sm;
	uint32_t* swp = sm + R / 2;
	uint32_t* sh = sm + R;
	const int q = blockIdx.y;
	const uint32_t p = t.p[q], pinv = t.p_inv[q];
	const uint32_t* wt = INV ? t.iw : t.w;
	const uint32_t* wtp = INV ? t.iwp : t.wp;
	const int stride = 4096 / R;
	for (int i = threadIdx.x; i < R / 2; i += NTT_THREADS)
	{
		sw[i] = wt[(size_t)q * 2048 + i * stride];
		swp[i] = wtp[(size_t)q * 2048 + i * stride];
	}
	const size_t blocks_per_row = C / COLS;
	const size_t row = blockIdx.x / blocks_per_row;
	const int lc = threadIdx.x % COLS, qq = threadIdx.x / COLS;
	const size_t l = (blockIdx.x % blocks_per_row) * COLS + lc;
	uint32_t* g = X + ((size_t)q * rows + row) * T + l;
	const uint32_t* TH = (INV ? t.iTH : t.TH) + ((size_t)q << (t.log_max - 11));
	const uint32_t* TL = (INV ? t.iTL : t.TL) + (size_t)q * 2048;
	const int tw_shift = t.log_max - logT; //omega_T = omega_Tmax^(2^tw_shift)
	const int npass = (logR + 2) / 3;

	uint32_t v[8];
	int idx[8], r, sb;
	pass_params(logR, 0, INV, r, sb);
	get_index_set(idx, r, qq, sb);
#pragma unroll
	for (int i = 0; i < 8; i++)
	{
		v[i] = g[(size_t)idx[i] * C];
	}
	if (INV) //untwist
	{
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			unsigned k1 = __brev((unsigned)idx[i]) >> (32 - logR);
			size_t E = ((size_t)l * k1) << tw_shift;
			v[i] = mul_mont(v[i], mul_mont(TH[E >> 11], TL[E & 2047], p, pinv), p, pinv);
		}
	}
	__syncthreads();
	for (int k = 0; k < npass; k++)
	{
		run_pass<INV>(v, r, qq, sb, R, sw, swp, p);
		if (k == npass - 1)
		{
			break;
		}
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			sh[COLS * (idx[i] + (idx[i] >> 3)) + lc] = v[i];
		}
		__syncthreads();
		pass_params(logR, k + 1, INV, r, sb);
		get_index_set(idx, r, qq, sb);
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			v[i] = sh[COLS * (idx[i] + (idx[i] >> 3)) + lc];
		}
		__syncthreads();
	}
	if (!INV) //twist
	{
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			unsigned k1 = __brev((unsigned)idx[i]) >> (32 - logR);
			size_t E = ((size_t)l * k1) << tw_shift;
			v[i] = mul_mont(v[i], mul_mont(TH[E >> 11], TL[E & 2047], p, pinv), p, pinv);
		}
	}
#pragma unroll
	for (int i = 0; i < 8; i++)
	{
		g[(size_t)idx[i] * C] = v[i];
	}
}

size_t row_smem(int logM)
{
	return (((size_t)1 << logM) + TILE + TILE / 8) * 4;
}

size_t col_smem(int logR)
{
	return (((size_t)1 << logR) + TILE + TILE / 8) * 4;
}

int split_logR(int logT)
{
	return (logT - 12 > 3) ? logT - 12 : 3;
}

unsigned int grid_size(size_t work, int threads)
{
	size_t blocks = (work + threads - 1) / threads;
	return (unsigned int)((blocks < 65536) ? ((blocks > 0) ? blocks : 1) : 65536);
}

template <bool INV>
void launch_row(const Gpu_ntt_tables& t, int primes, uint32_t* X, size_t nrows, int logM, cudaStream_t s)
{
	const size_t M = (size_t)1 << logM, rows_per_block = TILE / M;
	dim3 grid((unsigned)((nrows + rows_per_block - 1) / rows_per_block), (unsigned)primes);
	kRow<INV><<<grid, NTT_THREADS, row_smem(logM), s>>>(X, nrows, logM, t);
	CUDA_CHECK(cudaGetLastError());
}

template <bool INV>
void launch_col(const Gpu_ntt_tables& t, int primes, uint32_t* X, size_t rows, int logT, int logR, cudaStream_t s)
{
	const size_t C = (size_t)1 << (logT - logR), COLS = TILE >> logR;
	dim3 grid((unsigned)(rows * (C / COLS)), (unsigned)primes);
	kCol<INV><<<grid, NTT_THREADS, col_smem(logR), s>>>(X, rows, logT, logR, t);
	CUDA_CHECK(cudaGetLastError());
}

//The coefficients of b bits of the numbers of a level, modulo each prime, zero-padded to T.
__global__ void k_split(const uint64_t* __restrict__ limbs, size_t stride, const uint32_t* __restrict__ size, size_t rows, int b, int logT, uint32_t* __restrict__ X, Gpu_ntt_tables t)
{
	const int q = blockIdx.y;
	const uint32_t p = t.p[q];
	const size_t T = (size_t)1 << logT, total = rows * T;
	for (size_t e = (size_t)blockIdx.x * blockDim.x + threadIdx.x; e < total; e += (size_t)gridDim.x * blockDim.x)
	{
		const size_t row = e >> logT, i = e & (T - 1);
		const size_t n = size[row];
		const uint64_t* a = limbs + row * stride;
		const uint64_t pos = (uint64_t)i * b;
		const size_t li = (size_t)(pos >> 6);
		const unsigned int shift = (unsigned int)(pos & 63);
		uint64_t v = 0;
		if (li < n)
		{
			v = a[li] >> shift;
			if (shift != 0 && li + 1 < n)
			{
				v |= a[li + 1] << (64 - shift);
			}
			if (b < 64)
			{
				v &= (1ULL << b) - 1;
			}
		}
		uint32_t hi = (uint32_t)(v >> 32), lo = (uint32_t)v;
		while (lo >= p)
		{
			lo -= p;
		}
		X[((size_t)q * rows + row) * T + i] = add_mod(mul_shoup(hi, t.two32[q], t.two32_shoup[q], p), lo, p);
	}
}

__global__ void k_pointwise(uint32_t* __restrict__ X, const uint32_t* __restrict__ Y, size_t n, Gpu_ntt_tables t)
{
	const int q = blockIdx.y;
	const uint32_t p = t.p[q], pinv = t.p_inv[q];
	const size_t off = (size_t)q * n;
	for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < n; i += (size_t)gridDim.x * blockDim.x)
	{
		X[off + i] = mul_mont(X[off + i], Y[off + i], p, pinv);
	}
}

//MODE 0: Y[j] = X[2j] X[2j + 1]; MODE 1: Z[2j] = B[j] C[2j + 1], Z[2j + 1] = B[j] C[2j].
template <int MODE>
__global__ void k_pointwise_pairs(const uint32_t* __restrict__ B, const uint32_t* __restrict__ C, uint32_t* __restrict__ Z, size_t count, size_t T, Gpu_ntt_tables t)
{
	const int q = blockIdx.y;
	const uint32_t p = t.p[q], pinv = t.p_inv[q];
	const size_t total = (MODE == 0 ? count : 2 * count) * T;
	for (size_t e = (size_t)blockIdx.x * blockDim.x + threadIdx.x; e < total; e += (size_t)gridDim.x * blockDim.x)
	{
		const size_t row = e / T, pos = e % T;
		if (MODE == 0)
		{
			const uint32_t* x = B + ((size_t)q * 2 * count + 2 * row) * T + pos;
			Z[((size_t)q * count + row) * T + pos] = mul_mont(x[0], x[T], p, pinv);
		}
		else
		{
			const size_t j = row / 2, partner = row ^ 1;
			Z[((size_t)q * 2 * count + row) * T + pos] = mul_mont(B[((size_t)q * count + j) * T + pos], C[((size_t)q * 2 * count + partner) * T + pos], p, pinv);
		}
	}
}

//The coefficients of rows inverse transforms by the CRT, as 160-bit numbers in 3 limbs: C[(row T + k) 3 + 0..2].
__global__ void k_crt(const uint32_t* __restrict__ X, size_t rows, int logT, int P, uint64_t* __restrict__ C, Gpu_ntt_tables t)
{
	const size_t T = (size_t)1 << logT, total = rows * T;
	for (size_t e = (size_t)blockIdx.x * blockDim.x + threadIdx.x; e < total; e += (size_t)gridDim.x * blockDim.x)
	{
		const size_t row = e >> logT, k = e & (T - 1);
		uint32_t v[NTT_MAX_PRIMES];
		for (int q = 0; q < P; q++)
		{
			v[q] = mul_shoup(X[((size_t)q * rows + row) * T + k], t.scale[q][logT], t.scale_shoup[q][logT], t.p[q]);
		}
		for (int i = 1; i < P; i++)
		{
			const uint32_t pi = t.p[i];
			uint32_t x = v[i];
			for (int j = 0; j < i; j++)
			{
				uint32_t vj = (v[j] >= pi) ? v[j] - pi : v[j];
				x = mul_shoup(sub_mod(x, vj, pi), t.garner[j][i], t.garner_shoup[j][i], pi);
			}
			v[i] = x;
		}
		uint32_t words[6] = { v[P - 1], 0, 0, 0, 0, 0 };
		for (int i = P - 2; i >= 0; i--)
		{
			uint64_t carry = v[i];
			for (int w = 0; w < P; w++)
			{
				uint64_t x = (uint64_t)words[w] * t.p[i] + carry;
				words[w] = (uint32_t)x;
				carry = x >> 32;
			}
		}
		uint64_t* c = C + e * 3;
		c[0] = words[0] | ((uint64_t)words[1] << 32);
		c[1] = words[2] | ((uint64_t)words[3] << 32);
		c[2] = words[4];
	}
}

constexpr int COMBINE_THREADS = 256;
constexpr int COMBINE_ITEMS = 8;

//The sum of the parts of the coefficients c_k 2^(k b) that fall in limb j, as lo + 2^64 hi.
__device__ void limb_sum(const uint64_t* __restrict__ c, size_t T, int b, size_t j, uint64_t& lo, uint64_t& hi)
{
	lo = 0;
	hi = 0;
	const long long pos = 64 * (long long)j;
	long long k_lo = (pos >= 160) ? (pos - 160) / b + 1 : 0;
	long long k_hi = (pos + 63) / b;
	if (k_hi > (long long)T - 1)
	{
		k_hi = (long long)T - 1;
	}
	for (long long k = k_lo; k <= k_hi; k++)
	{
		const long long d = k * b - pos;
		const uint64_t* v = c + 3 * k;
		uint64_t chunk;
		if (d >= 0)
		{
			chunk = v[0] << d;
		}
		else
		{
			const long long e = -d;
			const int w = (int)(e >> 6), s = (int)(e & 63);
			chunk = (w < 3) ? v[w] >> s : 0;
			if (s != 0 && w + 1 < 3)
			{
				chunk |= v[w + 1] << (64 - s);
			}
		}
		lo += chunk;
		hi += (lo < chunk);
	}
}

//A carry function in 2 bits: bit c is the carry out for a carry in of c.
__device__ __forceinline__ uint32_t carry_apply(uint32_t f, uint32_t carry)
{
	return (f >> carry) & 1;
}

//f, then g.
__device__ __forceinline__ uint32_t carry_compose(uint32_t f, uint32_t g)
{
	return carry_apply(g, f & 1) | (carry_apply(g, (f >> 1) & 1) << 1);
}

constexpr uint32_t CARRY_IDENTITY = 2;

//Adds the inverse of the block's carries: each thread gives the carry function of its items, and gets the carry function of
//the items of the threads before it.
__device__ uint32_t block_exclusive_carry_scan(uint32_t f, uint32_t* scan)
{
	const int tid = threadIdx.x;
	scan[tid] = f;
	__syncthreads();
	for (int offset = 1; offset < COMBINE_THREADS; offset *= 2)
	{
		uint32_t v = (tid >= offset) ? carry_compose(scan[tid - offset], scan[tid]) : scan[tid];
		__syncthreads();
		scan[tid] = v;
		__syncthreads();
	}
	uint32_t prefix = (tid > 0) ? scan[tid - 1] : CARRY_IDENTITY;
	return prefix;
}

//The largest index below n of a nonzero limb, plus one (0 if all are zero), for the whole block.
__device__ size_t block_size(const uint64_t* o, size_t n, unsigned int* shared_max)
{
	if (threadIdx.x == 0)
	{
		*shared_max = 0;
	}
	__syncthreads();
	for (size_t j = threadIdx.x; j < n; j += blockDim.x)
	{
		if (o[j] != 0)
		{
			atomicMax(shared_max, (unsigned int)(j + 1));
		}
	}
	__syncthreads();
	size_t result = *shared_max;
	__syncthreads();
	return result;
}

//Adds a small number (3 limbs) at the bottom of n limbs, one thread. Returns the carry out.
__device__ uint64_t add_small(uint64_t* o, size_t n, const uint64_t add[3])
{
	uint64_t carry = 0;
	size_t j = 0;
	for (; j < n && (j < 3 || carry != 0); j++)
	{
		uint64_t a = (j < 3) ? add[j] : 0;
		uint64_t s = o[j] + a;
		uint64_t c = (s < a);
		s += carry;
		c |= (s < carry);
		o[j] = s;
		carry = c;
	}
	return carry;
}

//One block per row: the overlap-add of the coefficients into limbs, with carries resolved by a block scan, chunk by chunk.
//For a cyclic product, the part from bit L upwards is then added onto the bottom, with the end-around carry, and 2^L - 1 is 0.
__global__ void __launch_bounds__(COMBINE_THREADS) k_combine(const uint64_t* __restrict__ C, int logT, int b, int cyclic, uint64_t* __restrict__ out, size_t stride, uint32_t* __restrict__ out_size)
{
	__shared__ uint32_t scan[COMBINE_THREADS];
	__shared__ uint32_t chunk_carry;
	__shared__ unsigned int shared_max;
	__shared__ int all_ones;
	const size_t T = (size_t)1 << logT;
	const size_t row = blockIdx.x;
	const uint64_t* c = C + row * T * 3;
	uint64_t* o = out + row * stride;
	const size_t n_out = (((T - 1) * (size_t)b + 160) >> 6) + 2;
	const int tid = threadIdx.x;

	if (tid == 0)
	{
		chunk_carry = 0;
	}
	__syncthreads();
	for (size_t base = 0; base < n_out; base += COMBINE_THREADS * COMBINE_ITEMS)
	{
		const size_t j0 = base + (size_t)tid * COMBINE_ITEMS;
		uint64_t r[COMBINE_ITEMS];
		uint32_t f[COMBINE_ITEMS];
		uint32_t composed = CARRY_IDENTITY;
		uint64_t lo, prev_hi = 0;
		if (j0 > 0 && j0 - 1 < n_out)
		{
			limb_sum(c, T, b, j0 - 1, lo, prev_hi);
		}
#pragma unroll
		for (int it = 0; it < COMBINE_ITEMS; it++)
		{
			const size_t j = j0 + it;
			r[it] = 0;
			f[it] = CARRY_IDENTITY;
			if (j < n_out)
			{
				uint64_t hi;
				limb_sum(c, T, b, j, lo, hi);
				uint64_t s = lo + prev_hi;
				uint32_t g = (s < lo);
				r[it] = s;
				f[it] = g | ((g | (s == ~0ULL ? 1u : 0u)) << 1);
				prev_hi = hi;
			}
			composed = carry_compose(composed, f[it]);
		}
		const uint32_t carry_in = chunk_carry;
		uint32_t prefix = block_exclusive_carry_scan(composed, scan);
		uint32_t carry = carry_apply(prefix, carry_in);
#pragma unroll
		for (int it = 0; it < COMBINE_ITEMS; it++)
		{
			const size_t j = j0 + it;
			if (j < n_out)
			{
				o[j] = r[it] + carry;
				carry = carry_apply(f[it], carry);
			}
		}
		__syncthreads();
		if (tid == COMBINE_THREADS - 1)
		{
			chunk_carry = carry_apply(scan[COMBINE_THREADS - 1], carry_in);
		}
		__syncthreads();
	}
	for (size_t j = n_out + tid; j < stride; j += COMBINE_THREADS)
	{
		o[j] = 0;
	}
	__syncthreads();

	if (cyclic)
	{
		const uint64_t L = (uint64_t)T * b;
		const size_t L_limbs = (size_t)(L >> 6);
		const unsigned int L_bit = (unsigned int)(L & 63);
		const uint64_t L_mask = (1ULL << L_bit) - 1;
		if (tid == 0)
		{
			//The part from bit L upwards is less than 2^160.
			uint64_t high[3];
			for (int w = 0; w < 3; w++)
			{
				size_t j = L_limbs + w;
				uint64_t v = (j < n_out) ? o[j] >> L_bit : 0;
				if (L_bit != 0 && j + 1 < n_out)
				{
					v |= o[j + 1] << (64 - L_bit);
				}
				high[w] = v;
			}
			o[L_limbs] &= L_mask;
			for (size_t j = L_limbs + 1; j < n_out; j++)
			{
				o[j] = 0;
			}
			add_small(o, L_limbs + 1, high);
			if ((o[L_limbs] >> L_bit) & 1)
			{
				o[L_limbs] &= L_mask;
				uint64_t one[3] = { 1, 0, 0 };
				add_small(o, L_limbs + 1, one);
			}
			all_ones = (L_bit == 0 || o[L_limbs] == L_mask);
		}
		__syncthreads();
		for (size_t j = tid; j < L_limbs; j += COMBINE_THREADS)
		{
			if (o[j] != ~0ULL)
			{
				all_ones = 0;
			}
		}
		__syncthreads();
		if (all_ones)
		{
			for (size_t j = tid; j <= L_limbs; j += COMBINE_THREADS)
			{
				o[j] = 0;
			}
		}
		__syncthreads();
	}

	size_t size = block_size(o, n_out, &shared_max);
	if (tid == 0)
	{
		out_size[row] = (uint32_t)size;
	}
}

//One block per row: out = bits shift up to precision of in, minus one, masked, or 0 for y = 0.
__global__ void k_extract(const uint64_t* __restrict__ in, size_t in_stride, const uint32_t* __restrict__ shift_bits, const uint32_t* __restrict__ precision_bits,
	const uint32_t* __restrict__ y_size, int y_div, uint64_t* __restrict__ out, size_t out_stride, uint32_t* __restrict__ out_size)
{
	__shared__ unsigned int shared_min;
	__shared__ unsigned int shared_max;
	const size_t row = blockIdx.x;
	const int tid = threadIdx.x;
	const uint64_t* a = in + row * in_stride;
	uint64_t* o = out + row * out_stride;
	if (y_size[row / y_div] == 0)
	{
		if (tid == 0)
		{
			out_size[row] = 0;
		}
		return;
	}
	const uint64_t shift = shift_bits[row], precision = precision_bits[row];
	const uint64_t out_bits = precision - shift;
	const size_t out_limbs = (size_t)((out_bits + 63) >> 6);
	const size_t first = (size_t)(shift >> 6);
	const size_t source = (size_t)((precision + 63) >> 6) - first;
	const unsigned int s = (unsigned int)(shift & 63);
	for (size_t j = tid; j < source; j += blockDim.x)
	{
		uint64_t v = a[first + j] >> s;
		if (s != 0 && j + 1 < source)
		{
			v |= a[first + j + 1] << (64 - s);
		}
		o[j] = v;
	}
	if (tid == 0)
	{
		shared_min = 0xFFFFFFFFu;
	}
	__syncthreads();
	for (size_t j = tid; j < out_limbs; j += blockDim.x)
	{
		if (o[j] != 0)
		{
			atomicMin(&shared_min, (unsigned int)j);
		}
	}
	__syncthreads();
	const size_t z = shared_min; //the first nonzero limb: the borrow of the subtraction of one stops there
	for (size_t j = tid; j < out_limbs; j += blockDim.x)
	{
		if (j < z)
		{
			o[j] = ~0ULL;
		}
		else if (j == z)
		{
			o[j] -= 1;
		}
	}
	__syncthreads();
	if (tid == 0 && (out_bits & 63) != 0)
	{
		o[out_limbs - 1] &= (1ULL << (out_bits & 63)) - 1;
	}
	__syncthreads();
	size_t size = block_size(o, out_limbs, &shared_max);
	if (tid == 0)
	{
		out_size[row] = (uint32_t)size;
	}
}

//out = a b (na + nb limbs, normalized size returned).
__device__ uint32_t school_multiply(const uint64_t* a, uint32_t na, const uint64_t* b, uint32_t nb, uint64_t* out)
{
	if (na == 0 || nb == 0)
	{
		return 0;
	}
	for (uint32_t k = 0; k < na + nb; k++)
	{
		out[k] = 0;
	}
	for (uint32_t i = 0; i < na; i++)
	{
		const uint64_t ai = a[i];
		uint64_t carry = 0;
		for (uint32_t j = 0; j < nb; j++)
		{
			uint64_t lo = ai * b[j];
			uint64_t hi = __umul64hi(ai, b[j]);
			uint64_t o = out[i + j];
			lo += o;
			hi += (lo < o);
			lo += carry;
			hi += (lo < carry);
			out[i + j] = lo;
			carry = hi;
		}
		out[i + nb] = carry;
	}
	uint32_t n = na + nb;
	while (n > 0 && out[n - 1] == 0)
	{
		n--;
	}
	return n;
}

__global__ void k_school_products(const uint64_t* __restrict__ x, size_t x_stride, const uint32_t* __restrict__ x_size, size_t count, uint64_t* __restrict__ out, size_t out_stride, uint32_t* __restrict__ out_size)
{
	for (size_t j = (size_t)blockIdx.x * blockDim.x + threadIdx.x; j < count; j += (size_t)gridDim.x * blockDim.x)
	{
		out_size[j] = school_multiply(x + 2 * j * x_stride, x_size[2 * j], x + (2 * j + 1) * x_stride, x_size[2 * j + 1], out + j * out_stride);
	}
}

__global__ void k_school_middle(const uint64_t* __restrict__ y, size_t y_stride, const uint32_t* __restrict__ y_size, const uint64_t* __restrict__ x, size_t x_stride,
	const uint32_t* __restrict__ x_size, const uint32_t* __restrict__ shift_bits, const uint32_t* __restrict__ precision_bits, size_t rows,
	uint64_t* __restrict__ scratch, uint64_t* __restrict__ out, size_t out_stride, uint32_t* __restrict__ out_size)
{
	for (size_t r = (size_t)blockIdx.x * blockDim.x + threadIdx.x; r < rows; r += (size_t)gridDim.x * blockDim.x)
	{
		const size_t partner = r ^ 1;
		uint64_t* product = scratch + r * (y_stride + x_stride);
		uint32_t n = school_multiply(y + (r / 2) * y_stride, y_size[r / 2], x + partner * x_stride, x_size[partner], product);
		const uint64_t shift = shift_bits[r], precision = precision_bits[r];
		uint64_t* o = out + r * out_stride;
		//bits shift up to precision of the product, exactly
		const size_t first = (size_t)(shift >> 6);
		const unsigned int s = (unsigned int)(shift & 63);
		const uint64_t out_bits = precision - shift;
		const size_t out_limbs = (size_t)((out_bits + 63) >> 6);
		for (size_t j = 0; j < out_limbs; j++)
		{
			const size_t k = first + j;
			uint64_t v = (k < n) ? product[k] >> s : 0;
			if (s != 0 && k + 1 < n)
			{
				v |= product[k + 1] << (64 - s);
			}
			o[j] = v;
		}
		if ((out_bits & 63) != 0)
		{
			o[out_limbs - 1] &= (1ULL << (out_bits & 63)) - 1;
		}
		size_t size = out_limbs;
		while (size > 0 && o[size - 1] == 0)
		{
			size--;
		}
		out_size[r] = (uint32_t)size;
	}
}

__global__ void k_bits(const uint64_t* __restrict__ limbs, size_t stride, const uint32_t* __restrict__ size, size_t count, uint32_t* __restrict__ bits)
{
	for (size_t r = (size_t)blockIdx.x * blockDim.x + threadIdx.x; r < count; r += (size_t)gridDim.x * blockDim.x)
	{
		const uint32_t n = size[r];
		bits[r] = (n == 0) ? 1 : 64 * (n - 1) + (64 - __clzll((long long)limbs[r * stride + n - 1]));
	}
}

} //namespace

Gpu_ntt::Gpu_ntt()
{
	const Ntt_tables& host = ntt_tables();
	const size_t high_count = ((size_t)1 << NTT_MAX_LOG_LENGTH) / 2048;
	std::vector<uint32_t> all;
	//w, wp, iw, iwp, TL, iTL: [prime][2048]; TH, iTH: [prime][high_count]
	const size_t per_table_low = (size_t)NTT_MAX_PRIMES * 2048, per_table_high = (size_t)NTT_MAX_PRIMES * high_count;
	all.resize(6 * per_table_low + 2 * per_table_high);
	uint32_t* w = all.data();
	uint32_t* wp = w + per_table_low;
	uint32_t* iw = wp + per_table_low;
	uint32_t* iwp = iw + per_table_low;
	uint32_t* TL = iwp + per_table_low;
	uint32_t* iTL = TL + per_table_low;
	uint32_t* TH = iTL + per_table_low;
	uint32_t* iTH = TH + per_table_high;

	tables.log_max = NTT_MAX_LOG_LENGTH;
	for (int q = 0; q < NTT_MAX_PRIMES; q++)
	{
		const Ntt_prime& prime = host.prime[q];
		tables.p[q] = prime.p;
		tables.p_inv[q] = prime.p_inv;
		tables.two32[q] = prime.two32;
		tables.two32_shoup[q] = prime.two32_shoup;
		for (int lt = 0; lt <= NTT_MAX_LOG_LENGTH; lt++)
		{
			tables.scale[q][lt] = host.scale[q][lt];
			tables.scale_shoup[q][lt] = host.scale_shoup[q][lt];
		}
		for (int j = 0; j < NTT_MAX_PRIMES; j++)
		{
			tables.garner[j][q] = host.garner[j][q];
			tables.garner_shoup[j][q] = host.garner_shoup[j][q];
		}
		//The CPU tables hold omega_(2s)^j at s + j; omega_4096^t is at 2048 + t.
		for (int t = 0; t < 2048; t++)
		{
			w[q * 2048 + t] = prime.w[2048 + t];
			wp[q * 2048 + t] = prime.w_shoup[2048 + t];
			iw[q * 2048 + t] = prime.iw[2048 + t];
			iwp[q * 2048 + t] = prime.iw_shoup[2048 + t];
			TL[q * 2048 + t] = prime.twist_low[t];
			iTL[q * 2048 + t] = prime.itwist_low[t];
		}
		for (size_t h = 0; h < high_count; h++)
		{
			TH[q * high_count + h] = prime.twist_high[h];
			iTH[q * high_count + h] = prime.itwist_high[h];
		}
	}
	uint32_t* device;
	CUDA_CHECK(cudaMalloc(&device, all.size() * sizeof(uint32_t)));
	CUDA_CHECK(cudaMemcpy(device, all.data(), all.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));
	table_memory = device;
	tables.w = device;
	tables.wp = device + per_table_low;
	tables.iw = device + 2 * per_table_low;
	tables.iwp = device + 3 * per_table_low;
	tables.TL = device + 4 * per_table_low;
	tables.iTL = device + 5 * per_table_low;
	tables.TH = device + 6 * per_table_low;
	tables.iTH = device + 6 * per_table_low + per_table_high;

	scratch = NULL;
	scratch_limbs = 0;

	//The transforms of up to 4096 points need more than the default 48 KB of shared memory.
	CUDA_CHECK(cudaFuncSetAttribute(kRow<false>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)row_smem(12)));
	CUDA_CHECK(cudaFuncSetAttribute(kRow<true>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)row_smem(12)));
	CUDA_CHECK(cudaFuncSetAttribute(kCol<false>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)col_smem(12)));
	CUDA_CHECK(cudaFuncSetAttribute(kCol<true>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)col_smem(12)));
}

Gpu_ntt::~Gpu_ntt()
{
	cudaFree(table_memory);
	cudaFree(scratch);
}

void Gpu_ntt::Reserve_scratch(size_t limbs)
{
	if (limbs > scratch_limbs)
	{
		CUDA_CHECK(cudaDeviceSynchronize());
		cudaFree(scratch);
		CUDA_CHECK(cudaMalloc(&scratch, limbs * sizeof(uint64_t)));
		scratch_limbs = limbs;
	}
}

size_t Gpu_ntt::Combine_limbs(const Ntt_plan& plan)
{
	return (((plan.Length() - 1) * (size_t)plan.coefficient_bits + 160) >> 6) + 2;
}

void Gpu_ntt::Transform(const Ntt_plan& plan, uint32_t* X, size_t rows, cudaStream_t stream)
{
	const int logT = plan.log_length;
	if (logT <= 12)
	{
		launch_row<false>(tables, plan.primes, X, rows, logT, stream);
	}
	else
	{
		const int logR = split_logR(logT);
		launch_col<false>(tables, plan.primes, X, rows, logT, logR, stream);
		launch_row<false>(tables, plan.primes, X, rows << logR, logT - logR, stream);
	}
}

void Gpu_ntt::Inverse(const Ntt_plan& plan, uint32_t* X, size_t rows, cudaStream_t stream)
{
	const int logT = plan.log_length;
	if (logT <= 12)
	{
		launch_row<true>(tables, plan.primes, X, rows, logT, stream);
	}
	else
	{
		const int logR = split_logR(logT);
		launch_row<true>(tables, plan.primes, X, rows << logR, logT - logR, stream);
		launch_col<true>(tables, plan.primes, X, rows, logT, logR, stream);
	}
}

void Gpu_ntt::Forward(const Ntt_plan& plan, const Gpu_level& in, uint32_t* X, cudaStream_t stream)
{
	if (in.count == 0)
	{
		return;
	}
	dim3 grid(grid_size(in.count * plan.Length(), 256), (unsigned)plan.primes);
	k_split<<<grid, 256, 0, stream>>>(in.limbs, in.stride, in.size, in.count, plan.coefficient_bits, plan.log_length, X, tables);
	CUDA_CHECK(cudaGetLastError());
	Transform(plan, X, in.count, stream);
}

void Gpu_ntt::Pointwise(const Ntt_plan& plan, uint32_t* X, const uint32_t* Y, size_t rows, cudaStream_t stream)
{
	if (rows == 0)
	{
		return;
	}
	const size_t n = rows * plan.Length();
	k_pointwise<<<dim3(grid_size(n, 256), (unsigned)plan.primes), 256, 0, stream>>>(X, Y, n, tables);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::Pointwise_adjacent(const Ntt_plan& plan, const uint32_t* X, uint32_t* Y, size_t count, cudaStream_t stream)
{
	if (count == 0)
	{
		return;
	}
	k_pointwise_pairs<0><<<dim3(grid_size(count * plan.Length(), 256), (unsigned)plan.primes), 256, 0, stream>>>(X, NULL, Y, count, plan.Length(), tables);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::Pointwise_swap(const Ntt_plan& plan, const uint32_t* B, const uint32_t* C, uint32_t* Z, size_t count, cudaStream_t stream)
{
	if (count == 0)
	{
		return;
	}
	k_pointwise_pairs<1><<<dim3(grid_size(2 * count * plan.Length(), 256), (unsigned)plan.primes), 256, 0, stream>>>(B, C, Z, count, plan.Length(), tables);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::Combine(const Ntt_plan& plan, uint32_t* X, size_t rows, bool cyclic, Gpu_level& out, cudaStream_t stream)
{
	if (rows == 0)
	{
		return;
	}
	if (out.stride < Combine_limbs(plan))
	{
		fprintf(stderr, "Gpu_ntt::Combine: output stride %zu is less than %zu limbs\n", out.stride, Combine_limbs(plan));
		exit(3);
	}
	Inverse(plan, X, rows, stream);
	Reserve_scratch(rows * plan.Length() * 3);
	k_crt<<<grid_size(rows * plan.Length(), 256), 256, 0, stream>>>(X, rows, plan.log_length, plan.primes, scratch, tables);
	CUDA_CHECK(cudaGetLastError());
	k_combine<<<(unsigned)rows, COMBINE_THREADS, 0, stream>>>(scratch, plan.log_length, plan.coefficient_bits, cyclic ? 1 : 0, out.limbs, out.stride, out.size);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::Extract_middle(const Gpu_level& in, const uint32_t* shift, const uint32_t* precision, const uint32_t* y_size, int y_div, Gpu_level& out, cudaStream_t stream)
{
	if (in.count == 0)
	{
		return;
	}
	k_extract<<<(unsigned)in.count, 256, 0, stream>>>(in.limbs, in.stride, shift, precision, y_size, y_div, out.limbs, out.stride, out.size);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::School_products(const Gpu_level& x, Gpu_level& out, cudaStream_t stream)
{
	const size_t count = x.count / 2;
	if (count == 0)
	{
		return;
	}
	k_school_products<<<grid_size(count, 128), 128, 0, stream>>>(x.limbs, x.stride, x.size, count, out.limbs, out.stride, out.size);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::School_middle(const Gpu_level& y, const Gpu_level& x, const uint32_t* shift, const uint32_t* precision, Gpu_level& out, cudaStream_t stream)
{
	if (x.count == 0)
	{
		return;
	}
	Reserve_scratch(x.count * (y.stride + x.stride));
	k_school_middle<<<grid_size(x.count, 128), 128, 0, stream>>>(y.limbs, y.stride, y.size, x.limbs, x.stride, x.size, shift, precision, x.count, scratch, out.limbs, out.stride, out.size);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::Bits(const Gpu_level& in, uint32_t* bits, cudaStream_t stream)
{
	if (in.count == 0)
	{
		return;
	}
	k_bits<<<grid_size(in.count, 256), 256, 0, stream>>>(in.limbs, in.stride, in.size, in.count, bits);
	CUDA_CHECK(cudaGetLastError());
}
