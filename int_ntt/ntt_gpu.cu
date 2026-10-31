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

//The transform kernels take THREADS threads, of 8 elements each: a tile of 8 THREADS elements per block. 256 threads (tiles of 2048)
//where the rows or columns have at most 2048 points, which gives more and smaller blocks per multiprocessor than 512 (the kernels
//are limited by latency, and their registers allow about 1024 threads per multiprocessor either way); 512 for rows and columns of
//4096 points.

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

//Where the first pass of kRow or kCol takes its input. LOAD_X: from X, in place. Otherwise X is only written, and the input
//comes from where the pass before it would have written it, so that the data goes through memory once less:
//  LOAD_SPLIT: the coefficients of b bits of the numbers of a level (forward transforms, as k_split);
//  LOAD_PLAIN, LOAD_ADJACENT, LOAD_SWAP: the pointwise products of transforms (inverse transforms, as Gpu_products).
enum { LOAD_X, LOAD_SPLIT, LOAD_PLAIN, LOAD_ADJACENT, LOAD_SWAP };

struct Ntt_source
{
	const uint64_t* limbs; //LOAD_SPLIT: the level
	size_t stride;
	const uint32_t* size;
	int b;
	const uint32_t* B; //the products: Z[r] = B[r] C[r], Z[j] = C[2j] C[2j + 1], or Z[2j] = B[j] C[2j + 1] and Z[2j + 1] = B[j] C[2j]
	const uint32_t* C;
	int log_rows; //kRow: the rows of each transform, 2^log_rows (the R of the four-step method, or 1)
};

//Coefficient i of b bits of the number a of n limbs (0 above it) modulo prime q, as a residue.
__device__ __forceinline__ uint32_t split_value(const uint64_t* a, size_t n, uint64_t i, int b, int q, const Gpu_ntt_tables& t)
{
	const uint32_t p = t.p[q];
	const uint64_t pos = i * (uint64_t)b;
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
	return add_mod(mul_shoup(hi, t.two32[q], t.two32_shoup[q], p), lo, p);
}

//Element e of row row (of nrows, prime q) of the products of LOAD_PLAIN, LOAD_ADJACENT or LOAD_SWAP, rows of M elements.
template <int LOAD>
__device__ __forceinline__ uint32_t product_value(const Ntt_source& src, size_t nrows, size_t row, size_t M, int e, int q, uint32_t p, uint32_t pinv)
{
	const size_t sub = row & (((size_t)1 << src.log_rows) - 1), r = row >> src.log_rows;
	if (LOAD == LOAD_PLAIN)
	{
		const size_t at = ((size_t)q * nrows + row) * M + e;
		return mul_mont(src.B[at], src.C[at], p, pinv);
	}
	if (LOAD == LOAD_ADJACENT)
	{
		const size_t at = ((size_t)q * 2 * nrows + ((2 * r) << src.log_rows) + sub) * M + e;
		return mul_mont(src.C[at], src.C[at + (M << src.log_rows)], p, pinv);
	}
	const size_t b_at = ((size_t)q * (nrows / 2) + ((r >> 1) << src.log_rows) + sub) * M + e;
	const size_t c_at = ((size_t)q * nrows + ((r ^ 1) << src.log_rows) + sub) * M + e;
	return mul_mont(src.B[b_at], src.C[c_at], p, pinv);
}

//M-point transforms of contiguous rows: X[(q nrows + row) M + i], with the input of the first pass as LOAD.
template <bool INV, int LOAD, int THREADS>
__global__ void __launch_bounds__(THREADS) kRow(uint32_t* __restrict__ X, size_t nrows, int logM, Gpu_ntt_tables t, Ntt_source src)
{
	extern __shared__ uint32_t sm[];
	const int M = 1 << logM, tpr = M / 8, rows_per_block = 8 * THREADS / M;
	uint32_t* sw = sm;
	uint32_t* swp = sm + M / 2;
	uint32_t* sh = sm + M;
	const int q = blockIdx.y;
	const uint32_t p = t.p[q];
	const uint32_t* wt = INV ? t.iw : t.w;
	const uint32_t* wtp = INV ? t.iwp : t.wp;
	const int stride = 4096 / M;
	for (int i = threadIdx.x; i < M / 2; i += THREADS)
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
		if (LOAD == LOAD_X)
		{
#pragma unroll
			for (int i = 0; i < 8; i++)
			{
				v[i] = g[idx[i]];
			}
		}
		else if (LOAD == LOAD_SPLIT)
		{
			const uint64_t* a = src.limbs + row * src.stride;
			const size_t n = src.size[row];
#pragma unroll
			for (int i = 0; i < 8; i++)
			{
				v[i] = split_value(a, n, (uint64_t)idx[i], src.b, q, t);
			}
		}
		else
		{
			const uint32_t pinv = t.p_inv[q];
#pragma unroll
			for (int i = 0; i < 8; i++)
			{
				v[i] = product_value<LOAD>(src, nrows, row, M, idx[i], q, p, pinv);
			}
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
//inverse pass). COLS = 4096 / R columns per block, R / 8 threads per column. With SPLIT (forward only), the input is split from
//the numbers of the level of src, as k_split.
template <bool INV, bool SPLIT, int THREADS>
__global__ void __launch_bounds__(THREADS) kCol(uint32_t* __restrict__ X, size_t rows, int logT, int logR, Gpu_ntt_tables t, Ntt_source src)
{
	extern __shared__ uint32_t sm[];
	const int R = 1 << logR, COLS = 8 * THREADS / R, logC = logT - logR;
	const size_t C = (size_t)1 << logC, T = (size_t)1 << logT;
	uint32_t* sw = sm;
	uint32_t* swp = sm + R / 2;
	uint32_t* sh = sm + R;
	const int q = blockIdx.y;
	const uint32_t p = t.p[q], pinv = t.p_inv[q];
	const uint32_t* wt = INV ? t.iw : t.w;
	const uint32_t* wtp = INV ? t.iwp : t.wp;
	const int stride = 4096 / R;
	for (int i = threadIdx.x; i < R / 2; i += THREADS)
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
	if (SPLIT)
	{
		const uint64_t* a = src.limbs + row * src.stride;
		const size_t n = src.size[row];
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			v[i] = split_value(a, n, (uint64_t)idx[i] * C + l, src.b, q, t);
		}
	}
	else
	{
#pragma unroll
		for (int i = 0; i < 8; i++)
		{
			v[i] = g[(size_t)idx[i] * C];
		}
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

size_t row_smem(int logM, int threads)
{
	return (((size_t)1 << logM) + 9 * (size_t)threads) * 4;
}

size_t col_smem(int logR, int threads)
{
	return (((size_t)1 << logR) + 9 * (size_t)threads) * 4;
}

//Transforms of up to 4096 points are one row (one kernel). Longer ones are split into rows of 2048 points and columns, of 256
//threads (ntt_core splits them into rows of 4096, but the result is the same, in bit-reversed order), or for the longest (more
//than 2^22), into rows of 4096.
bool single_row(int logT)
{
	return logT <= 12;
}

int split_logR(int logT)
{
	if (logT > 22)
	{
		return logT - 12;
	}
	return (logT - 11 > 3) ? logT - 11 : 3;
}

unsigned int grid_size(size_t work, int threads)
{
	size_t blocks = (work + threads - 1) / threads;
	return (unsigned int)((blocks < 65536) ? ((blocks > 0) ? blocks : 1) : 65536);
}

template <bool INV, int LOAD>
void launch_row(const Gpu_ntt_tables& t, int primes, uint32_t* X, size_t nrows, int logM, const Ntt_source& src, cudaStream_t s)
{
	const int threads = (logM <= 11) ? 256 : 512;
	const size_t M = (size_t)1 << logM, rows_per_block = 8 * (size_t)threads / M;
	dim3 grid((unsigned)((nrows + rows_per_block - 1) / rows_per_block), (unsigned)primes);
	if (threads == 256)
	{
		kRow<INV, LOAD, 256><<<grid, 256, row_smem(logM, 256), s>>>(X, nrows, logM, t, src);
	}
	else
	{
		kRow<INV, LOAD, 512><<<grid, 512, row_smem(logM, 512), s>>>(X, nrows, logM, t, src);
	}
	CUDA_CHECK(cudaGetLastError());
}

template <bool INV, bool SPLIT>
void launch_col(const Gpu_ntt_tables& t, int primes, uint32_t* X, size_t rows, int logT, int logR, const Ntt_source& src, cudaStream_t s)
{
	const int threads = (logR <= 11) ? 256 : 512;
	const size_t C = (size_t)1 << (logT - logR), COLS = (8 * (size_t)threads) >> logR;
	dim3 grid((unsigned)(rows * (C / COLS)), (unsigned)primes);
	if (threads == 256)
	{
		kCol<INV, SPLIT, 256><<<grid, 256, col_smem(logR, 256), s>>>(X, rows, logT, logR, t, src);
	}
	else
	{
		kCol<INV, SPLIT, 512><<<grid, 512, col_smem(logR, 512), s>>>(X, rows, logT, logR, t, src);
	}
	CUDA_CHECK(cudaGetLastError());
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

//Coefficient k of row row of the inverse transforms, by the CRT, as a 160-bit number in 5 words.
//P primes, a template parameter so that the loops are unrolled and the residues stay in registers.
template <int P>
__device__ __forceinline__ void crt_coefficient(const uint32_t* __restrict__ X, size_t rows, int logT, size_t row, size_t k, const Gpu_ntt_tables& t, uint32_t* words)
{
	const size_t T = (size_t)1 << logT;
	uint32_t v[P];
#pragma unroll
	for (int q = 0; q < P; q++)
	{
		v[q] = mul_shoup(X[((size_t)q * rows + row) * T + k], t.scale[q][logT], t.scale_shoup[q][logT], t.p[q]);
	}
#pragma unroll
	for (int i = 1; i < P; i++)
	{
		const uint32_t pi = t.p[i];
		uint32_t x = v[i];
#pragma unroll
		for (int j = 0; j < i; j++)
		{
			uint32_t vj = (v[j] >= pi) ? v[j] - pi : v[j];
			x = mul_shoup(sub_mod(x, vj, pi), t.garner[j][i], t.garner_shoup[j][i], pi);
		}
		v[i] = x;
	}
	uint32_t w[5] = { v[P - 1], 0, 0, 0, 0 };
#pragma unroll
	for (int i = P - 2; i >= 0; i--)
	{
		uint64_t carry = v[i];
#pragma unroll
		for (int k2 = 0; k2 < P; k2++)
		{
			uint64_t x = (uint64_t)w[k2] * t.p[i] + carry;
			w[k2] = (uint32_t)x;
			carry = x >> 32;
		}
	}
#pragma unroll
	for (int k2 = 0; k2 < 5; k2++)
	{
		words[k2] = w[k2];
	}
}

//The coefficients of rows inverse transforms by the CRT, as 160-bit numbers in 3 limbs: C[(row T + k) 3 + 0..2].
template <int P>
__global__ void k_crt(const uint32_t* __restrict__ X, size_t rows, int logT, uint64_t* __restrict__ C, Gpu_ntt_tables t)
{
	const size_t T = (size_t)1 << logT, total = rows * T;
	for (size_t e = (size_t)blockIdx.x * blockDim.x + threadIdx.x; e < total; e += (size_t)gridDim.x * blockDim.x)
	{
		uint32_t words[5];
		crt_coefficient<P>(X, rows, logT, e >> logT, e & (T - 1), t, words);
		uint64_t* c = C + e * 3;
		c[0] = words[0] | ((uint64_t)words[1] << 32);
		c[1] = words[2] | ((uint64_t)words[3] << 32);
		c[2] = words[4];
	}
}

constexpr int COMBINE_THREADS = 256;
constexpr int COMBINE_ITEMS = 4;

//The coefficients from which limb j is summed: from k_lo to k_hi.
__device__ __forceinline__ void limb_coefficients(size_t T, int b, size_t j, long long& k_lo, long long& k_hi)
{
	const long long pos = 64 * (long long)j;
	k_lo = (pos >= 160) ? (pos - 160) / b + 1 : 0;
	k_hi = (pos + 63) / b;
	if (k_hi > (long long)T - 1)
	{
		k_hi = (long long)T - 1;
	}
}

//The sum of the parts of the coefficients c_k 2^(k b) that fall in limb j, as lo + 2^64 hi. coefficient(k, v) gives c_k in 3 limbs.
template <class Coefficient>
__device__ __forceinline__ void limb_sum(Coefficient coefficient, size_t T, int b, size_t j, uint64_t& lo, uint64_t& hi)
{
	lo = 0;
	hi = 0;
	const long long pos = 64 * (long long)j;
	long long k_lo, k_hi;
	limb_coefficients(T, b, j, k_lo, k_hi);
	for (long long k = k_lo; k <= k_hi; k++)
	{
		const long long d = k * b - pos;
		uint64_t v[3];
		coefficient(k, v);
		uint64_t chunk;
		if (d >= 0)
		{
			chunk = v[0] << d;
		}
		else
		{
			const long long e = -d;
			const int w = (int)(e >> 6), s = (int)(e & 63);
			const uint64_t low = (w == 0) ? v[0] : ((w == 1) ? v[1] : ((w == 2) ? v[2] : 0));
			const uint64_t high = (w == 0) ? v[1] : ((w == 1) ? v[2] : 0);
			chunk = low >> s;
			if (s != 0)
			{
				chunk |= high << (64 - s);
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

//The carry function of a chunk (chunk_f) also has bit 2 set if a limb of the chunk generates a carry (or a borrow). Without it,
//and without a carry into the chunk, the carries change none of its limbs. The compositions only look at bits 0 and 1.
constexpr uint32_t CHUNK_GENERATES = 4;

//The largest x of the warp, in all its lanes.
__device__ __forceinline__ unsigned int warp_max(unsigned int x)
{
	for (int offset = 16; offset > 0; offset /= 2)
	{
		const unsigned int y = __shfl_xor_sync(0xFFFFFFFFu, x, offset);
		x = (y > x) ? y : x;
	}
	return x;
}

//The exclusive scan of the carry functions of a block: each thread gives the carry function of its items, and gets the carry
//function of the items of the threads before it. scan then holds the inclusive scan.
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
	//From the top down, blockDim limbs at a time: the top nonzero limb is usually in the first of them.
	for (size_t top = n; top > 0;)
	{
		const size_t width = (top < blockDim.x) ? top : blockDim.x;
		if (threadIdx.x < width && __ldcg(o + top - 1 - threadIdx.x) != 0)
		{
			atomicMax(shared_max, (unsigned int)(top - threadIdx.x));
		}
		__syncthreads();
		const unsigned int found = *shared_max;
		__syncthreads();
		if (found != 0)
		{
			break;
		}
		top -= width;
	}
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

//The combine of rows of CRT coefficients into integers, in four kernels, so that a row of millions of limbs is spread over many
//blocks: the overlap-add of each chunk of COMBINE_THREADS * COMBINE_ITEMS limbs, with the carry function of each limb and of the
//chunk; the carries into each chunk; the carries applied; and per row, the wraparound of a cyclic product and the size.
constexpr size_t COMBINE_CHUNK = (size_t)COMBINE_THREADS * COMBINE_ITEMS;

__host__ __device__ __forceinline__ size_t combine_limbs(int logT, int b)
{
	return (((((size_t)1 << logT) - 1) * (size_t)b + 160) >> 6) + 2;
}

//The coefficients a chunk of limbs needs at most, for coefficients of b bits.
__host__ __device__ __forceinline__ size_t chunk_coefficients(int b)
{
	return (64 * (COMBINE_CHUNK + 1) + 223) / (size_t)b + 2;
}

//The most coefficients that k_combine_limbs<true> keeps in shared memory, 5 words each, if the GPU allows that much shared
//memory per block (Gpu_ntt::combine_coefficients).
constexpr size_t COMBINE_SHARED_COEFFICIENTS = 4096;

//With FUSED, the coefficients are computed here by the CRT from the inverse transforms X, the ones of the chunk into shared
//memory, instead of being read from C (k_crt). That saves writing them, and reading each of them several times.
template <bool FUSED, int P>
__global__ void __launch_bounds__(COMBINE_THREADS) k_combine_limbs(const uint64_t* __restrict__ C, const uint32_t* __restrict__ X, size_t rows, Gpu_ntt_tables t,
	int logT, int b, size_t chunks, uint64_t* __restrict__ out, size_t stride, uint8_t* __restrict__ limb_f, uint8_t* __restrict__ chunk_f, unsigned int* __restrict__ row_info)
{
	__shared__ uint32_t scan[COMBINE_THREADS];
	extern __shared__ uint32_t coefficient_words[];
	const size_t T = (size_t)1 << logT;
	const size_t row = blockIdx.x / chunks, chunk = blockIdx.x % chunks;
	const uint64_t* c = C + row * T * 3;
	uint64_t* o = out + row * stride;
	const size_t n_out = combine_limbs(logT, b);
	uint8_t* f_row = limb_f + row * n_out;
	const size_t j0 = chunk * COMBINE_CHUNK + (size_t)threadIdx.x * COMBINE_ITEMS;

	long long first = 0;
	if (FUSED)
	{
		//The coefficients of the limbs of the chunk, and of the limb before it.
		const size_t j_first = chunk * COMBINE_CHUNK;
		size_t j_last = j_first + COMBINE_CHUNK - 1;
		j_last = (j_last < n_out - 1) ? j_last : n_out - 1;
		long long lo_k, hi_k, unused;
		limb_coefficients(T, b, (j_first > 0) ? j_first - 1 : 0, lo_k, unused);
		limb_coefficients(T, b, j_last, unused, hi_k);
		first = lo_k;
		for (long long k = lo_k + threadIdx.x; k <= hi_k; k += COMBINE_THREADS)
		{
			crt_coefficient<P>(X, rows, logT, row, (size_t)k, t, coefficient_words + 5 * (k - lo_k));
		}
		__syncthreads();
	}
	auto coefficient = [&](long long k, uint64_t v[3])
	{
		if (FUSED)
		{
			const uint32_t* w = coefficient_words + 5 * (k - first);
			v[0] = w[0] | ((uint64_t)w[1] << 32);
			v[1] = w[2] | ((uint64_t)w[3] << 32);
			v[2] = w[4];
		}
		else
		{
			v[0] = c[3 * k];
			v[1] = c[3 * k + 1];
			v[2] = c[3 * k + 2];
		}
	};

	uint32_t composed = CARRY_IDENTITY, generates = 0;
	unsigned int top = 0;
	uint64_t lo, prev_hi = 0;
	if (j0 > 0 && j0 - 1 < n_out)
	{
		limb_sum(coefficient, T, b, j0 - 1, lo, prev_hi);
	}
	for (int it = 0; it < COMBINE_ITEMS; it++)
	{
		const size_t j = j0 + it;
		if (j < n_out)
		{
			uint64_t hi;
			limb_sum(coefficient, T, b, j, lo, hi);
			uint64_t s = lo + prev_hi;
			uint32_t g = (s < lo);
			uint32_t f = g | ((g | (s == ~0ULL ? 1u : 0u)) << 1);
			o[j] = s;
			f_row[j] = (uint8_t)f;
			composed = carry_compose(composed, f);
			generates |= g;
			top = (s != 0) ? (unsigned int)(j + 1) : top;
			prev_hi = hi;
		}
	}
	block_exclusive_carry_scan(composed, scan);
	const int any_generates = __syncthreads_or(generates);
	if (threadIdx.x == COMBINE_THREADS - 1)
	{
		chunk_f[row * chunks + chunk] = (uint8_t)(scan[COMBINE_THREADS - 1] | (any_generates ? CHUNK_GENERATES : 0));
	}
	//The top nonzero limb of the row before the carries, for its size (row_info[2 row + 1]).
	top = warp_max(top);
	if ((threadIdx.x & 31) == 0 && top != 0)
	{
		atomicMax(row_info + 2 * row + 1, top);
	}
}

//o += 2^(64 from), modulo 2^(64 (last + 1)), by the whole block: the limbs of all ones from limb from on become 0, and the first
//other one up to limb last gains 1. shared is a word of shared memory.
__device__ void block_increment(uint64_t* o, size_t from, size_t last, unsigned int* shared)
{
	if (threadIdx.x == 0)
	{
		*shared = ~0u;
	}
	__syncthreads();
	for (size_t base = from; base <= last; base += blockDim.x)
	{
		const size_t j = base + threadIdx.x;
		if (j <= last && __ldcg(o + j) != ~0ULL)
		{
			atomicMin(shared, (unsigned int)j);
		}
		__syncthreads();
		const unsigned int found = *shared;
		__syncthreads();
		if (found != ~0u)
		{
			break;
		}
	}
	const size_t z = *shared;
	for (size_t j = from + threadIdx.x; j < z && j <= last; j += blockDim.x)
	{
		o[j] = 0;
	}
	if (threadIdx.x == 0 && z <= last)
	{
		o[z] = __ldcg(o + z) + 1;
	}
	__syncthreads();
}

//The finish of a combined row, by one block: the limbs above the product are zeroed. For a cyclic product, the part from bit L
//upwards is added onto the bottom, with the end-around carry, and 2^L - 1 is 0. Then the size: top, the top nonzero limb plus one
//found by the blocks of the row, for an exact product; for a cyclic one, from there or from the limb at bit L down, whichever is
//lower (without a part from bit L upwards, nothing above top changes). The limbs written by other blocks of the same kernel are
//read through L2 (__ldcg).
__device__ void combine_finish(int logT, int b, int cyclic, uint64_t* __restrict__ o, size_t stride, uint32_t* __restrict__ out_size, unsigned int top)
{
	__shared__ unsigned int shared_max;
	__shared__ int all_ones, carry_on;
	const size_t T = (size_t)1 << logT;
	const size_t n_out = combine_limbs(logT, b);
	const int tid = threadIdx.x;

	for (size_t j = n_out + tid; j < stride; j += blockDim.x)
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
				uint64_t v = (j < n_out) ? __ldcg(o + j) >> L_bit : 0;
				if (L_bit != 0 && j + 1 < n_out)
				{
					v |= __ldcg(o + j + 1) << (64 - L_bit);
				}
				high[w] = v;
			}
			o[L_limbs] = __ldcg(o + L_limbs) & L_mask;
			for (size_t j = L_limbs + 1; j < n_out; j++)
			{
				o[j] = 0;
			}
			uint64_t carry = 0;
			for (size_t j = 0; j <= L_limbs && j < 3; j++)
			{
				const uint64_t a = high[j];
				uint64_t s = __ldcg(o + j) + a;
				uint64_t c = (s < a);
				s += carry;
				c |= (s < carry);
				o[j] = s;
				carry = c;
			}
			carry_on = (carry != 0 && L_limbs >= 3);
		}
		__syncthreads();
		//The carry can run through many limbs of all ones (the root works with complements), so the whole block takes it on.
		if (carry_on)
		{
			block_increment(o, 3, L_limbs, &shared_max);
		}
		if (tid == 0)
		{
			carry_on = (__ldcg(o + L_limbs) >> L_bit) & 1;
			if (carry_on)
			{
				o[L_limbs] &= L_mask;
			}
		}
		__syncthreads();
		//The end-around carry.
		if (carry_on)
		{
			block_increment(o, 0, L_limbs, &shared_max);
		}
		if (tid == 0)
		{
			all_ones = (L_bit == 0 || __ldcg(o + L_limbs) == L_mask);
		}
		__syncthreads();
		//2^L - 1, all ones: checked blockDim limbs at a time, stopping at the first that is not.
		for (size_t base = 0; base < L_limbs; base += blockDim.x)
		{
			const size_t j = base + tid;
			if (j < L_limbs && __ldcg(o + j) != ~0ULL)
			{
				all_ones = 0;
			}
			__syncthreads();
			if (!all_ones)
			{
				break;
			}
		}
		if (all_ones)
		{
			for (size_t j = tid; j <= L_limbs; j += blockDim.x)
			{
				o[j] = 0;
			}
		}
		__syncthreads();
	}

	if (!cyclic)
	{
		if (tid == 0)
		{
			out_size[0] = top;
		}
		return;
	}
	const size_t L_limbs = (size_t)(((uint64_t)T * b) >> 6);
	const size_t from = (L_limbs + 1 < (size_t)top) ? L_limbs + 1 : (size_t)top;
	size_t size = block_size(o, (from < n_out) ? from : n_out, &shared_max);
	if (tid == 0)
	{
		out_size[0] = (uint32_t)size;
	}
}

//The carries of rows of n limbs, from the carry function of each limb (limb_f) and of each chunk of COMBINE_CHUNK limbs (chunk_f),
//with carry_in into each row. Each block works out the carry into its chunk from the chunks before it, and applies the carries,
//if there are any in the chunk (they are rare: a limb sum must overflow). row_info has two words per row: the blocks done, and the
//top nonzero limb plus one (from the combine, and raised here by a carry), both 0 again after the row. The last block of a row to
//finish then writes the carry out of the row, or for a subtraction the borrow, if carry_out is not NULL, and does the finish of a
//Combine if finish is set.
__global__ void __launch_bounds__(COMBINE_THREADS) k_carry(size_t n, size_t chunks, uint64_t* __restrict__ out, size_t stride, const uint8_t* __restrict__ limb_f,
	const uint8_t* __restrict__ chunk_f, uint32_t carry_in, uint32_t* __restrict__ carry_out, unsigned int* __restrict__ row_info, int finish, int logT, int b, int cyclic,
	uint32_t* __restrict__ out_size)
{
	__shared__ uint32_t scan[COMBINE_THREADS];
	__shared__ int last;
	const int tid = threadIdx.x;
	const size_t row = blockIdx.x / chunks, chunk = blockIdx.x % chunks;
	uint64_t* o = out + row * stride;
	const uint8_t* f_row = limb_f + row * n;
	const uint8_t* chunks_row = chunk_f + row * chunks;

	//The carry into this chunk: the chunks before it, a range of them for each thread.
	const size_t per = (chunk + COMBINE_THREADS - 1) / COMBINE_THREADS;
	uint32_t before = CARRY_IDENTITY;
	for (size_t k = tid * per; k < (tid + 1) * per && k < chunk; k++)
	{
		before = carry_compose(before, chunks_row[k]);
	}
	block_exclusive_carry_scan(before, scan);
	const uint32_t chunk_in = carry_apply(scan[COMBINE_THREADS - 1], carry_in);
	__syncthreads();

	if (chunk_in != 0 || (chunks_row[chunk] & CHUNK_GENERATES))
	{
		const size_t j0 = chunk * COMBINE_CHUNK + (size_t)tid * COMBINE_ITEMS;
		uint32_t f[COMBINE_ITEMS];
		uint32_t composed = CARRY_IDENTITY;
#pragma unroll
		for (int it = 0; it < COMBINE_ITEMS; it++)
		{
			const size_t j = j0 + it;
			f[it] = (j < n) ? f_row[j] : CARRY_IDENTITY;
			composed = carry_compose(composed, f[it]);
		}
		uint32_t prefix = block_exclusive_carry_scan(composed, scan);
		uint32_t carry = carry_apply(prefix, chunk_in);
		unsigned int top = 0;
#pragma unroll
		for (int it = 0; it < COMBINE_ITEMS; it++)
		{
			const size_t j = j0 + it;
			if (j < n)
			{
				if (carry != 0)
				{
					const uint64_t v = o[j] + carry;
					o[j] = v;
					top = (v != 0) ? (unsigned int)(j + 1) : top;
				}
				carry = carry_apply(f[it], carry);
			}
		}
		//A carry can reach above the top nonzero limb.
		top = warp_max(top);
		if ((tid & 31) == 0 && top != 0)
		{
			atomicMax(row_info + 2 * row + 1, top);
		}
	}

	if (carry_out == NULL && !finish)
	{
		return;
	}
	__threadfence();
	__syncthreads();
	if (tid == 0)
	{
		last = (atomicAdd(row_info + 2 * row, 1u) == (unsigned int)(chunks - 1));
	}
	__syncthreads();
	if (!last)
	{
		return;
	}
	__shared__ unsigned int row_top;
	if (tid == 0)
	{
		row_info[2 * row] = 0;
		row_top = atomicExch(row_info + 2 * row + 1, 0u);
	}
	__syncthreads();
	if (carry_out != NULL)
	{
		//The carry out of the row: all its chunks.
		const size_t all = (chunks + COMBINE_THREADS - 1) / COMBINE_THREADS;
		uint32_t total = CARRY_IDENTITY;
		for (size_t k = tid * all; k < (tid + 1) * all && k < chunks; k++)
		{
			total = carry_compose(total, chunks_row[k]);
		}
		__syncthreads();
		block_exclusive_carry_scan(total, scan);
		if (tid == 0)
		{
			carry_out[row] = carry_apply(scan[COMBINE_THREADS - 1], carry_in) ^ carry_in;
		}
	}
	if (finish)
	{
		combine_finish(logT, b, cyclic, o, stride, out_size + row, row_top);
	}
}

//Extract_middle in two kernels, so that both the few long rows at the top of a tree and the many short ones at the bottom are
//spread over the GPU. out = bits shift up to precision of in, minus one, masked, or 0 for y = 0.
constexpr int EXTRACT_THREADS = 128;
constexpr int EXTRACT_ITEMS = 4;
constexpr size_t EXTRACT_CHUNK = (size_t)EXTRACT_THREADS * EXTRACT_ITEMS;

//The bits, a block per chunk of EXTRACT_CHUNK limbs of a row, and the lowest nonzero limb of each row (first, which is ~0 before).
__global__ void __launch_bounds__(EXTRACT_THREADS) k_extract_copy(const uint64_t* __restrict__ in, size_t in_stride, const uint32_t* __restrict__ shift_bits,
	const uint32_t* __restrict__ precision_bits, const uint32_t* __restrict__ y_size, int y_div, size_t chunks, uint64_t* __restrict__ out, size_t out_stride,
	unsigned int* __restrict__ first)
{
	const size_t row = blockIdx.x / chunks, chunk = blockIdx.x % chunks;
	if (y_size[row / y_div] == 0)
	{
		return;
	}
	const uint64_t shift = shift_bits[row], precision = precision_bits[row];
	const size_t out_limbs = (size_t)((precision - shift + 63) >> 6);
	const size_t first_limb = (size_t)(shift >> 6);
	const size_t source = (size_t)((precision + 63) >> 6) - first_limb;
	const unsigned int s = (unsigned int)(shift & 63);
	const uint64_t* a = in + row * in_stride + first_limb;
	uint64_t* o = out + row * out_stride;
	unsigned int lowest = ~0u;
#pragma unroll
	for (int it = 0; it < EXTRACT_ITEMS; it++)
	{
		const size_t j = chunk * EXTRACT_CHUNK + (size_t)it * EXTRACT_THREADS + threadIdx.x;
		if (j < source)
		{
			uint64_t v = a[j] >> s;
			if (s != 0 && j + 1 < source)
			{
				v |= a[j + 1] << (64 - s);
			}
			o[j] = v;
			if (v != 0 && j < out_limbs && j < lowest)
			{
				lowest = (unsigned int)j;
			}
		}
	}
	for (int offset = 16; offset > 0; offset /= 2)
	{
		const unsigned int other = __shfl_xor_sync(0xFFFFFFFFu, lowest, offset);
		lowest = (other < lowest) ? other : lowest;
	}
	if ((threadIdx.x & 31) == 0 && lowest != ~0u)
	{
		atomicMin(first + row, lowest);
	}
}

//Then a warp per row: minus one (the borrow stops at the lowest nonzero limb), the mask of the top limb, and the size. first is
//~0 again afterwards.
__global__ void k_extract_finish(const uint32_t* __restrict__ shift_bits, const uint32_t* __restrict__ precision_bits, const uint32_t* __restrict__ y_size, int y_div,
	size_t rows, uint64_t* __restrict__ out, size_t out_stride, uint32_t* __restrict__ out_size, unsigned int* __restrict__ first)
{
	const int lane = threadIdx.x & 31;
	const size_t row = ((size_t)blockIdx.x * blockDim.x + threadIdx.x) / 32;
	if (row >= rows)
	{
		return;
	}
	if (y_size[row / y_div] == 0)
	{
		if (lane == 0)
		{
			out_size[row] = 0;
		}
		return;
	}
	const uint64_t out_bits = precision_bits[row] - shift_bits[row];
	const size_t out_limbs = (size_t)((out_bits + 63) >> 6);
	uint64_t* o = out + row * out_stride;
	const size_t z = first[row];
	__syncwarp();
	if (lane == 0)
	{
		first[row] = ~0u;
	}
	for (size_t j = lane; j < out_limbs && j <= z; j += 32)
	{
		if (j < z)
		{
			o[j] = ~0ULL;
		}
		else
		{
			o[j] -= 1;
		}
	}
	__syncwarp();
	if (lane == 0 && (out_bits & 63) != 0)
	{
		o[out_limbs - 1] &= (1ULL << (out_bits & 63)) - 1;
	}
	__syncwarp();
	//The size, from the top down, 32 limbs at a time: lane i looks at limb top - 1 - i.
	size_t size = 0;
	for (size_t top = out_limbs; top > 0; top = (top > 32) ? top - 32 : 0)
	{
		const bool nonzero = (size_t)lane < top && o[top - 1 - lane] != 0;
		const unsigned int found = __ballot_sync(0xFFFFFFFFu, nonzero);
		if (found != 0)
		{
			size = top - (size_t)(__ffs(found) - 1);
			break;
		}
	}
	if (lane == 0)
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

//out = a + b or a - b (as a + ~b + 1) over n limbs, chunk by chunk as k_combine_limbs: the sums without carries in, the carry
//function of each limb, and of each chunk. b is b_limbs limbs at limb b_offset, or if b is NULL, the limb b_value at b_offset.
__global__ void __launch_bounds__(COMBINE_THREADS) k_add_limbs(const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs, size_t b_offset, uint64_t b_value,
	int subtract, size_t n, uint64_t* out, uint8_t* __restrict__ limb_f, uint8_t* __restrict__ chunk_f)
{
	__shared__ uint32_t scan[COMBINE_THREADS];
	const size_t j0 = (size_t)blockIdx.x * COMBINE_CHUNK + (size_t)threadIdx.x * COMBINE_ITEMS;
	uint32_t composed = CARRY_IDENTITY, generates = 0;
	for (int it = 0; it < COMBINE_ITEMS; it++)
	{
		const size_t j = j0 + it;
		if (j < n)
		{
			const uint64_t x = (a != NULL && j < a_limbs) ? a[j] : 0;
			uint64_t y;
			if (b != NULL)
			{
				y = (j >= b_offset && j - b_offset < b_limbs) ? b[j - b_offset] : 0;
			}
			else
			{
				y = (j == b_offset) ? b_value : 0;
			}
			y = subtract ? ~y : y;
			const uint64_t s = x + y;
			const uint32_t g = (s < x);
			const uint32_t f = g | ((g | (s == ~0ULL ? 1u : 0u)) << 1);
			out[j] = s;
			limb_f[j] = (uint8_t)f;
			composed = carry_compose(composed, f);
			generates |= g;
		}
	}
	block_exclusive_carry_scan(composed, scan);
	const int any_generates = __syncthreads_or(generates);
	if (threadIdx.x == COMBINE_THREADS - 1)
	{
		chunk_f[blockIdx.x] = (uint8_t)(scan[COMBINE_THREADS - 1] | (any_generates ? CHUNK_GENERATES : 0));
	}
}

//out = floor(in / 2^position) modulo 2^count (count 0: no limit), of in or of its complement (~in, with ones above in_limbs).
//A negative position shifts to the left.
__global__ void k_shift(const uint64_t* __restrict__ in, size_t in_limbs, long long position, uint64_t count, int complement, uint64_t* __restrict__ out, size_t out_limbs)
{
	const uint64_t flip = complement ? ~0ULL : 0;
	for (size_t j = (size_t)blockIdx.x * blockDim.x + threadIdx.x; j < out_limbs; j += (size_t)gridDim.x * blockDim.x)
	{
		const long long bit = 64 * (long long)j + position;
		const long long w = (bit >= 0) ? bit / 64 : -((-bit + 63) / 64);
		const unsigned int s = (unsigned int)(bit - 64 * w);
		uint64_t lo = (w >= 0 && w < (long long)in_limbs) ? in[w] ^ flip : (w >= 0 ? flip : 0);
		uint64_t v = lo >> s;
		if (s != 0)
		{
			const long long w1 = w + 1;
			uint64_t hi = (w1 >= 0 && w1 < (long long)in_limbs) ? in[w1] ^ flip : (w1 >= 0 ? flip : 0);
			v |= hi << (64 - s);
		}
		if (count != 0)
		{
			const uint64_t bottom = 64 * (uint64_t)j;
			if (bottom >= count)
			{
				v = 0;
			}
			else if (count - bottom < 64)
			{
				v &= (1ULL << (count - bottom)) - 1;
			}
		}
		out[j] = v;
	}
}

#define SIZE_CHUNK 4096 //limbs per block of k_size

//*limbs = the number of limbs of in without leading zeros, and *bits its size in bits (1 for 0), if bits is not NULL. Each block
//takes SIZE_CHUNK limbs, the first block the top ones, and skips them if a block above has found a nonzero limb already: the
//numbers often have many zero limbs on top, which one block would go through one step at a time. found[0] is the largest size
//found, and found[1] counts the blocks done: the last one writes the results and sets both to 0 for the next call.
__global__ void k_size(const uint64_t* __restrict__ in, size_t n, uint32_t* __restrict__ limbs, uint32_t* __restrict__ bits, unsigned int* found)
{
	__shared__ unsigned int shared_max;
	__shared__ bool skip, last;
	const size_t top = n - (size_t)blockIdx.x * SIZE_CHUNK;
	const size_t bottom = (top > SIZE_CHUNK) ? top - SIZE_CHUNK : 0;
	if (threadIdx.x == 0)
	{
		skip = *(volatile unsigned int*)found > top;
	}
	__syncthreads();
	if (!skip)
	{
		const size_t size = block_size(in + bottom, top - bottom, &shared_max);
		if (threadIdx.x == 0 && size != 0)
		{
			atomicMax(found, (unsigned int)(bottom + size));
		}
	}
	if (threadIdx.x == 0)
	{
		__threadfence();
		last = (atomicAdd(found + 1, 1u) == gridDim.x - 1);
	}
	__syncthreads();
	if (last && threadIdx.x == 0)
	{
		const size_t size = atomicExch(found, 0u);
		found[1] = 0;
		*limbs = (uint32_t)size;
		if (bits != NULL)
		{
			*bits = (size == 0) ? 1 : (uint32_t)(64 * (size - 1) + (64 - __clzll((long long)__ldcg(in + size - 1))));
		}
	}
}

__global__ void k_set_word(uint32_t* p, uint32_t value)
{
	*p = value;
}

//Small numbers in one thread, as arrays of limbs with a size.
__device__ __forceinline__ int limbs_normalize(const uint64_t* a, int n)
{
	while (n > 0 && a[n - 1] == 0)
	{
		n--;
	}
	return n;
}

//out = bits shift up to shift + count of a (n limbs), in out_limbs = ceil(count / 64) limbs. Returns the size without leading zeros.
__device__ int limbs_bits(const uint64_t* a, int n, uint64_t shift, uint64_t count, uint64_t* out)
{
	const int out_limbs = (int)((count + 63) >> 6);
	const int first = (int)(shift >> 6);
	const unsigned int s = (unsigned int)(shift & 63);
	for (int j = 0; j < out_limbs; j++)
	{
		const int k = first + j;
		uint64_t v = (k < n) ? a[k] >> s : 0;
		if (s != 0 && k + 1 < n)
		{
			v |= a[k + 1] << (64 - s);
		}
		out[j] = v;
	}
	if ((count & 63) != 0 && out_limbs > 0)
	{
		out[out_limbs - 1] &= (1ULL << (count & 63)) - 1;
	}
	return limbs_normalize(out, out_limbs);
}

//a = a >> s in place, for a nonzero shift; returns the new size.
__device__ int limbs_shift_right(uint64_t* a, int n, unsigned int s)
{
	const int w = (int)(s >> 6);
	const unsigned int b = s & 63;
	for (int j = 0; j + w < n; j++)
	{
		uint64_t v = a[j + w] >> b;
		if (b != 0 && j + w + 1 < n)
		{
			v |= a[j + w + 1] << (64 - b);
		}
		a[j] = v;
	}
	for (int j = (n - w > 0) ? n - w : 0; j < n; j++)
	{
		a[j] = 0;
	}
	return limbs_normalize(a, n);
}

//a = a << s, into n limbs (the result must fit); returns the new size.
__device__ int limbs_shift_left(uint64_t* a, int na, int n, unsigned int s)
{
	const int w = (int)(s >> 6);
	const unsigned int b = s & 63;
	for (int j = n - 1; j >= 0; j--)
	{
		const int k = j - w;
		uint64_t v = (k >= 0 && k < na) ? a[k] << b : 0;
		if (b != 0 && k - 1 >= 0 && k - 1 < na)
		{
			v |= a[k - 1] >> (64 - b);
		}
		a[j] = v;
	}
	return limbs_normalize(a, n);
}

__device__ unsigned int limbs_ctz(const uint64_t* a, int n)
{
	for (int j = 0; j < n; j++)
	{
		if (a[j] != 0)
		{
			return 64 * j + (unsigned int)(__ffsll((long long)a[j]) - 1);
		}
	}
	return 0;
}

__device__ int limbs_compare(const uint64_t* a, int na, const uint64_t* b, int nb)
{
	if (na != nb)
	{
		return (na > nb) ? 1 : -1;
	}
	for (int j = na - 1; j >= 0; j--)
	{
		if (a[j] != b[j])
		{
			return (a[j] > b[j]) ? 1 : -1;
		}
	}
	return 0;
}

//a = a - b, for a >= b; returns the new size.
__device__ int limbs_subtract(uint64_t* a, int na, const uint64_t* b, int nb)
{
	uint64_t borrow = 0;
	for (int j = 0; j < na; j++)
	{
		const uint64_t y = (j < nb) ? b[j] : 0;
		const uint64_t d = a[j] - y;
		const uint64_t borrow2 = (a[j] < y) | (d < borrow);
		a[j] = d - borrow;
		borrow = borrow2;
	}
	return limbs_normalize(a, na);
}

//q = n / g for odd g dividing n, by Hensel division from the bottom limb; n is destroyed. Returns the size of q.
__device__ int limbs_divexact_odd(uint64_t* n, int nn, const uint64_t* g, int ng, uint64_t* q)
{
	uint64_t inv = g[0]; //the inverse of g[0] modulo 2^64, by Newton's iteration from 3 bits
	for (int i = 0; i < 5; i++)
	{
		inv *= 2 - g[0] * inv;
	}
	for (int j = 0; j < nn; j++)
	{
		const uint64_t qj = n[j] * inv;
		q[j] = qj;
		//n -= qj g B^j
		uint64_t carry = 0, borrow = 0;
		for (int k = 0; j + k < nn; k++)
		{
			uint64_t lo = 0, hi = 0;
			if (k < ng)
			{
				lo = qj * g[k];
				hi = __umul64hi(qj, g[k]);
			}
			lo += carry;
			hi += (lo < carry);
			carry = hi;
			const uint64_t x = n[j + k];
			const uint64_t d = x - lo;
			const uint64_t borrow2 = (x < lo) | (d < borrow);
			n[j + k] = d - borrow;
			borrow = borrow2;
			if (k >= ng && carry == 0 && borrow == 0)
			{
				break;
			}
		}
	}
	return limbs_normalize(q, nn);
}

constexpr int LEAF_MAX_SQUARINGS = 2;

//The cofactor of each leaf N of the scaled remainder tree, one thread per leaf: from the fraction y of its parent (scaled by
//2^(parent bits + guard)) and its sibling x, r = P mod N = round(((y x) mod 2^(parent bits + guard)) >> (parent bits - N bits)
//times N / 2^(N bits + guard)), as Cpu_tree_backend::Leaf_cofactors, and then N / gcd(N, r^(2^squarings)). K is the most limbs
//of a leaf, and a parent fraction has at most 2K + 1.
template <int K>
__global__ void __launch_bounds__(128) k_leaf_cofactors(const uint64_t* __restrict__ y, size_t y_stride, const uint32_t* __restrict__ y_size, const uint32_t* __restrict__ y_bits,
	const uint64_t* __restrict__ x, size_t x_stride, const uint32_t* __restrict__ x_size, const uint32_t* __restrict__ x_bits, size_t count, uint32_t guard, int squarings,
	uint64_t* __restrict__ out, size_t out_stride, uint32_t* __restrict__ out_size)
{
	constexpr int V = K << LEAF_MAX_SQUARINGS; //r^(2^squarings)
	for (size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x; i < count; i += (size_t)gridDim.x * blockDim.x)
	{
		uint64_t product[3 * K + 1], q[K + 2], a[V], b[V], c[K + 1];
		const size_t parent = i / 2, sibling = i ^ 1;
		const uint32_t nb = x_bits[i], pb = y_bits[parent];
		const uint64_t* N = x + i * x_stride;
		const int n_size = (int)x_size[i];

		//r = P mod N
		int size = (int)school_multiply(y + parent * y_stride, y_size[parent], x + sibling * x_stride, x_size[sibling], product);
		size = limbs_bits(product, size, pb - nb, nb + guard, q);
		size = (int)school_multiply(q, size, N, n_size, product);
		size = limbs_bits(product, size, nb, guard + 64 * (K + 1), q);
		uint64_t carry = (1ULL << guard) - 1;
		for (int j = 0; j < K + 2 && carry != 0; j++)
		{
			q[j] += carry;
			carry = (q[j] < carry);
		}
		size = limbs_bits(q, K + 2, guard, 64 * K, a);

		//v = r^(2^squarings), in b
		for (int s = 0; s < squarings; s++)
		{
			size = (int)school_multiply(a, size, a, size, b);
			for (int j = 0; j < size; j++)
			{
				a[j] = b[j];
			}
		}
		for (int j = 0; j < size; j++)
		{
			b[j] = a[j];
		}
		int vb = size;

		//g = gcd(N, v) by the binary algorithm, in a
		for (int j = 0; j < n_size; j++)
		{
			a[j] = N[j];
		}
		int va = n_size;
		if (vb != 0 && va != 0)
		{
			const unsigned int za = limbs_ctz(a, va), zb = limbs_ctz(b, vb);
			const unsigned int z = (za < zb) ? za : zb;
			if (za != 0)
			{
				va = limbs_shift_right(a, va, za);
			}
			if (zb != 0)
			{
				vb = limbs_shift_right(b, vb, zb);
			}
			uint64_t *u = a, *v = b;
			int nu = va, nv = vb;
			for (;;)
			{
				//u and v odd
				const int order = limbs_compare(u, nu, v, nv);
				if (order == 0)
				{
					break;
				}
				if (order > 0)
				{
					uint64_t* t = u;
					u = v;
					v = t;
					int tn = nu;
					nu = nv;
					nv = tn;
				}
				nv = limbs_subtract(v, nv, u, nu);
				nv = limbs_shift_right(v, nv, limbs_ctz(v, nv));
			}
			if (u != a)
			{
				for (int j = 0; j < nu; j++)
				{
					a[j] = u[j];
				}
			}
			va = nu;
			if (z != 0)
			{
				for (int j = va; j < K + 1 && j < V; j++)
				{
					a[j] = 0;
				}
				va = limbs_shift_left(a, va, K, z);
			}
		}
		//else gcd(N, 0) = N, already in a

		//c = N / g: the common power of 2 removed, then Hensel division by the odd part.
		for (int j = 0; j < n_size; j++)
		{
			b[j] = N[j];
		}
		int nn = n_size;
		const unsigned int zg = limbs_ctz(a, va);
		if (zg != 0)
		{
			va = limbs_shift_right(a, va, zg);
			nn = limbs_shift_right(b, nn, zg);
		}
		const int c_size = limbs_divexact_odd(b, nn, a, va, c);
		uint64_t* o = out + i * out_stride;
		for (int j = 0; j < K; j++)
		{
			o[j] = (j < c_size) ? c[j] : 0;
		}
		out_size[i] = (uint32_t)c_size;
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
	carry_scratch = NULL;
	carry_bytes = 0;
	row_done = NULL;
	row_done_count = 0;
	extract_first = NULL;
	extract_first_count = 0;
	CUDA_CHECK(cudaMalloc(&size_found, 2 * sizeof(unsigned int)));
	CUDA_CHECK(cudaMemset(size_found, 0, 2 * sizeof(unsigned int)));

	//The transforms of up to 4096 points need more than the default 48 KB of shared memory.
	CUDA_CHECK(cudaFuncSetAttribute(kRow<false, LOAD_X, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)row_smem(12, 512)));
	CUDA_CHECK(cudaFuncSetAttribute(kRow<false, LOAD_SPLIT, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)row_smem(12, 512)));
	CUDA_CHECK(cudaFuncSetAttribute(kRow<true, LOAD_X, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)row_smem(12, 512)));
	CUDA_CHECK(cudaFuncSetAttribute(kRow<true, LOAD_PLAIN, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)row_smem(12, 512)));
	CUDA_CHECK(cudaFuncSetAttribute(kRow<true, LOAD_ADJACENT, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)row_smem(12, 512)));
	CUDA_CHECK(cudaFuncSetAttribute(kRow<true, LOAD_SWAP, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)row_smem(12, 512)));
	CUDA_CHECK(cudaFuncSetAttribute(kCol<false, false, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)col_smem(12, 512)));
	CUDA_CHECK(cudaFuncSetAttribute(kCol<false, true, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)col_smem(12, 512)));
	CUDA_CHECK(cudaFuncSetAttribute(kCol<true, false, 512>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)col_smem(12, 512)));
	//The shared coefficients of the combine: as many as the GPU allows per block next to the rest of the block's shared memory
	//(64 KB on compute capability 7.5, 99 KB on 8.6 and 8.9), up to COMBINE_SHARED_COEFFICIENTS.
	int device_index = 0, optin = 48 * 1024;
	CUDA_CHECK(cudaGetDevice(&device_index));
	CUDA_CHECK(cudaDeviceGetAttribute(&optin, cudaDevAttrMaxSharedMemoryPerBlockOptin, device_index));
	combine_coefficients = ((size_t)optin - 4096) / (5 * sizeof(uint32_t));
	combine_coefficients = (combine_coefficients < COMBINE_SHARED_COEFFICIENTS) ? combine_coefficients : COMBINE_SHARED_COEFFICIENTS;
	const int combine_shared = (int)(5 * sizeof(uint32_t) * combine_coefficients);
	CUDA_CHECK(cudaFuncSetAttribute(k_combine_limbs<true, 3>, cudaFuncAttributeMaxDynamicSharedMemorySize, combine_shared));
	CUDA_CHECK(cudaFuncSetAttribute(k_combine_limbs<true, 4>, cudaFuncAttributeMaxDynamicSharedMemorySize, combine_shared));
	CUDA_CHECK(cudaFuncSetAttribute(k_combine_limbs<true, 5>, cudaFuncAttributeMaxDynamicSharedMemorySize, combine_shared));
}

Gpu_ntt::~Gpu_ntt()
{
	cudaFree(table_memory);
	cudaFree(scratch);
	cudaFree(carry_scratch);
	cudaFree(row_done);
	cudaFree(size_found);
	cudaFree(extract_first);
}

void Gpu_ntt::Reserve_scratch(size_t limbs)
{
	if (limbs > scratch_limbs)
	{
		//With some room to spare: the sizes vary a little, and cudaFree waits for the whole GPU.
		CUDA_CHECK(cudaDeviceSynchronize());
		cudaFree(scratch);
		scratch_limbs = limbs + limbs / 8;
		CUDA_CHECK(cudaMalloc(&scratch, scratch_limbs * sizeof(uint64_t)));
	}
}

void Gpu_ntt::Reserve_carry(size_t bytes)
{
	if (bytes > carry_bytes)
	{
		CUDA_CHECK(cudaDeviceSynchronize());
		cudaFree(carry_scratch);
		carry_bytes = bytes + bytes / 8;
		CUDA_CHECK(cudaMalloc(&carry_scratch, carry_bytes));
	}
}

size_t Gpu_ntt::Combine_limbs(const Ntt_plan& plan)
{
	return combine_limbs(plan.log_length, plan.coefficient_bits);
}

void Gpu_ntt::Transform(const Ntt_plan& plan, uint32_t* X, size_t rows, cudaStream_t stream)
{
	const int logT = plan.log_length;
	const Ntt_source none = {};
	if (single_row(logT))
	{
		launch_row<false, LOAD_X>(tables, plan.primes, X, rows, logT, none, stream);
	}
	else
	{
		const int logR = split_logR(logT);
		launch_col<false, false>(tables, plan.primes, X, rows, logT, logR, none, stream);
		launch_row<false, LOAD_X>(tables, plan.primes, X, rows << logR, logT - logR, none, stream);
	}
}

void Gpu_ntt::Inverse(const Ntt_plan& plan, uint32_t* X, size_t rows, cudaStream_t stream)
{
	const int logT = plan.log_length;
	const Ntt_source none = {};
	if (single_row(logT))
	{
		launch_row<true, LOAD_X>(tables, plan.primes, X, rows, logT, none, stream);
	}
	else
	{
		const int logR = split_logR(logT);
		launch_row<true, LOAD_X>(tables, plan.primes, X, rows << logR, logT - logR, none, stream);
		launch_col<true, false>(tables, plan.primes, X, rows, logT, logR, none, stream);
	}
}

//The split of the numbers into coefficients is done by the first pass of the transform.
void Gpu_ntt::Forward(const Ntt_plan& plan, const Gpu_level& in, uint32_t* X, cudaStream_t stream)
{
	if (in.count == 0)
	{
		return;
	}
	const int logT = plan.log_length;
	Ntt_source src = {};
	src.limbs = in.limbs;
	src.stride = in.stride;
	src.size = in.size;
	src.b = plan.coefficient_bits;
	if (single_row(logT))
	{
		launch_row<false, LOAD_SPLIT>(tables, plan.primes, X, in.count, logT, src, stream);
	}
	else
	{
		const int logR = split_logR(logT);
		const Ntt_source none = {};
		launch_col<false, true>(tables, plan.primes, X, in.count, logT, logR, src, stream);
		launch_row<false, LOAD_X>(tables, plan.primes, X, in.count << logR, logT - logR, none, stream);
	}
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
	Combine_inverse(plan, X, rows, cyclic, out, stream);
}

//The products are done by the first pass of the inverse transform, which reads them from B and C instead of Z.
void Gpu_ntt::Combine_products(const Ntt_plan& plan, Gpu_products kind, const uint32_t* B, const uint32_t* C, uint32_t* Z, size_t rows, bool cyclic, Gpu_level& out,
	cudaStream_t stream)
{
	if (rows == 0)
	{
		return;
	}
	if (out.stride < Combine_limbs(plan))
	{
		fprintf(stderr, "Gpu_ntt::Combine_products: output stride %zu is less than %zu limbs\n", out.stride, Combine_limbs(plan));
		exit(3);
	}
	const int logT = plan.log_length;
	const int logR = single_row(logT) ? 0 : split_logR(logT);
	Ntt_source src = {};
	src.B = B;
	src.C = C;
	src.log_rows = logR;
	const size_t nrows = rows << logR;
	const int logM = logT - logR;
	switch (kind)
	{
	case PRODUCTS_PLAIN:
		launch_row<true, LOAD_PLAIN>(tables, plan.primes, Z, nrows, logM, src, stream);
		break;
	case PRODUCTS_ADJACENT:
		launch_row<true, LOAD_ADJACENT>(tables, plan.primes, Z, nrows, logM, src, stream);
		break;
	default:
		launch_row<true, LOAD_SWAP>(tables, plan.primes, Z, nrows, logM, src, stream);
		break;
	}
	if (!single_row(logT))
	{
		const Ntt_source none = {};
		launch_col<true, false>(tables, plan.primes, Z, rows, logT, logR, none, stream);
	}
	Combine_inverse(plan, Z, rows, cyclic, out, stream);
}

//The CRT and the carries of Combine, after the inverse transforms.
void Gpu_ntt::Combine_inverse(const Ntt_plan& plan, uint32_t* X, size_t rows, bool cyclic, Gpu_level& out, cudaStream_t stream)
{
	const size_t n_out = Combine_limbs(plan), chunks = (n_out + COMBINE_CHUNK - 1) / COMBINE_CHUNK;
	Reserve_carry(rows * (n_out + chunks));
	uint8_t* limb_f = carry_scratch;
	uint8_t* chunk_f = carry_scratch + rows * n_out;
	const int logT = plan.log_length, b = plan.coefficient_bits;
	Reserve_rows(rows);
	if (chunk_coefficients(b) <= combine_coefficients)
	{
		const size_t coefficients = (chunk_coefficients(b) < plan.Length()) ? chunk_coefficients(b) : plan.Length();
		const size_t shared = 5 * sizeof(uint32_t) * coefficients;
		const unsigned grid = (unsigned)(rows * chunks);
		switch (plan.primes)
		{
		case 3:
			k_combine_limbs<true, 3><<<grid, COMBINE_THREADS, shared, stream>>>(NULL, X, rows, tables, logT, b, chunks, out.limbs, out.stride, limb_f, chunk_f, row_done);
			break;
		case 4:
			k_combine_limbs<true, 4><<<grid, COMBINE_THREADS, shared, stream>>>(NULL, X, rows, tables, logT, b, chunks, out.limbs, out.stride, limb_f, chunk_f, row_done);
			break;
		default:
			k_combine_limbs<true, 5><<<grid, COMBINE_THREADS, shared, stream>>>(NULL, X, rows, tables, logT, b, chunks, out.limbs, out.stride, limb_f, chunk_f, row_done);
			break;
		}
	}
	else
	{
		Reserve_scratch(rows * plan.Length() * 3);
		const unsigned crt_grid = grid_size(rows * plan.Length(), 256);
		switch (plan.primes)
		{
		case 3:
			k_crt<3><<<crt_grid, 256, 0, stream>>>(X, rows, logT, scratch, tables);
			break;
		case 4:
			k_crt<4><<<crt_grid, 256, 0, stream>>>(X, rows, logT, scratch, tables);
			break;
		default:
			k_crt<5><<<crt_grid, 256, 0, stream>>>(X, rows, logT, scratch, tables);
			break;
		}
		CUDA_CHECK(cudaGetLastError());
		k_combine_limbs<false, 3><<<(unsigned)(rows * chunks), COMBINE_THREADS, 0, stream>>>(scratch, X, rows, tables, logT, b, chunks, out.limbs, out.stride, limb_f, chunk_f,
			row_done);
	}
	CUDA_CHECK(cudaGetLastError());
	k_carry<<<(unsigned)(rows * chunks), COMBINE_THREADS, 0, stream>>>(n_out, chunks, out.limbs, out.stride, limb_f, chunk_f, 0, NULL, row_done, 1, logT, b, cyclic ? 1 : 0, out.size);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::Extract_middle(const Gpu_level& in, const uint32_t* shift, const uint32_t* precision, const uint32_t* y_size, int y_div, Gpu_level& out, cudaStream_t stream)
{
	if (in.count == 0)
	{
		return;
	}
	if (in.count > extract_first_count)
	{
		CUDA_CHECK(cudaDeviceSynchronize());
		cudaFree(extract_first);
		extract_first_count = in.count + in.count / 8;
		CUDA_CHECK(cudaMalloc(&extract_first, extract_first_count * sizeof(unsigned int)));
		CUDA_CHECK(cudaMemset(extract_first, 0xFF, extract_first_count * sizeof(unsigned int)));
	}
	const size_t chunks = (out.stride + EXTRACT_CHUNK - 1) / EXTRACT_CHUNK;
	k_extract_copy<<<(unsigned)(in.count * chunks), EXTRACT_THREADS, 0, stream>>>(in.limbs, in.stride, shift, precision, y_size, y_div, chunks, out.limbs, out.stride,
		extract_first);
	CUDA_CHECK(cudaGetLastError());
	k_extract_finish<<<(unsigned)((in.count * 32 + 255) / 256), 256, 0, stream>>>(shift, precision, y_size, y_div, in.count, out.limbs, out.stride, out.size, extract_first);
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

void Gpu_ntt::Shift(uint64_t* out, size_t out_limbs, const uint64_t* in, size_t in_limbs, long long position, uint64_t count, bool complement, cudaStream_t stream)
{
	if (out_limbs == 0)
	{
		return;
	}
	k_shift<<<grid_size(out_limbs, 256), 256, 0, stream>>>(in, in_limbs, position, count, complement ? 1 : 0, out, out_limbs);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::Add(uint64_t* out, size_t n, const uint64_t* a, size_t a_limbs, const uint64_t* b, size_t b_limbs, size_t b_offset, uint64_t b_value, bool subtract,
	uint32_t* carry, cudaStream_t stream)
{
	if (n == 0)
	{
		return;
	}
	const size_t chunks = (n + COMBINE_CHUNK - 1) / COMBINE_CHUNK;
	Reserve_carry(n + chunks);
	uint8_t* limb_f = carry_scratch;
	uint8_t* chunk_f = carry_scratch + n;
	k_add_limbs<<<(unsigned)chunks, COMBINE_THREADS, 0, stream>>>(a, a_limbs, b, b_limbs, b_offset, b_value, subtract ? 1 : 0, n, out, limb_f, chunk_f);
	CUDA_CHECK(cudaGetLastError());
	Reserve_rows(1);
	k_carry<<<(unsigned)chunks, COMBINE_THREADS, 0, stream>>>(n, chunks, out, 0, limb_f, chunk_f, subtract ? 1 : 0, carry, row_done, 0, 0, 0, 0, NULL);
	CUDA_CHECK(cudaGetLastError());
}

void Gpu_ntt::Size(const uint64_t* in, size_t n, uint32_t* limbs, uint32_t* bits, cudaStream_t stream)
{
	const size_t blocks = (n + SIZE_CHUNK - 1) / SIZE_CHUNK;
	k_size<<<(unsigned)((blocks > 0) ? blocks : 1), 256, 0, stream>>>(in, n, limbs, bits, size_found);
	CUDA_CHECK(cudaGetLastError());
}

//The counters of the blocks done in each row, for k_carry, which leaves them at 0.
void Gpu_ntt::Reserve_rows(size_t rows)
{
	if (rows > row_done_count)
	{
		CUDA_CHECK(cudaDeviceSynchronize());
		cudaFree(row_done);
		row_done_count = rows + rows / 8;
		CUDA_CHECK(cudaMalloc(&row_done, 2 * row_done_count * sizeof(unsigned int)));
		CUDA_CHECK(cudaMemset(row_done, 0, 2 * row_done_count * sizeof(unsigned int)));
	}
}

void Gpu_ntt::Set_word(uint32_t* p, uint32_t value, cudaStream_t stream)
{
	k_set_word<<<1, 1, 0, stream>>>(p, value);
	CUDA_CHECK(cudaGetLastError());
}

int Gpu_ntt::Leaf_kernel_limbs(int max_limbs)
{
	static const int K[] = { 2, 3, 4, 6, 8 };
	for (int k : K)
	{
		if (max_limbs <= k)
		{
			return k;
		}
	}
	return 0;
}

bool Gpu_ntt::Leaf_cofactors(const Gpu_level& y, const uint32_t* y_bits, const Gpu_level& x, const uint32_t* x_bits, int max_limbs, uint32_t guard, int squarings, Gpu_level& out,
	cudaStream_t stream)
{
	const int K = Leaf_kernel_limbs(max_limbs);
	if (K == 0 || squarings < 0 || squarings > LEAF_MAX_SQUARINGS || guard == 0 || guard >= 64 || out.stride < (size_t)K)
	{
		return false;
	}
	if (x.count == 0)
	{
		return true;
	}
	const unsigned int grid = grid_size(x.count, 128);
#define LEAF_CASE(k) case k: k_leaf_cofactors<k><<<grid, 128, 0, stream>>>(y.limbs, y.stride, y.size, y_bits, x.limbs, x.stride, x.size, x_bits, x.count, guard, squarings, out.limbs, out.stride, out.size); break;
	switch (K)
	{
		LEAF_CASE(2)
		LEAF_CASE(3)
		LEAF_CASE(4)
		LEAF_CASE(6)
		LEAF_CASE(8)
	}
#undef LEAF_CASE
	CUDA_CHECK(cudaGetLastError());
	return true;
}
