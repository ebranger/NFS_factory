//The transforms of ntt_core.h on the CPU. They follow the kernels of gpuecm2's ntt.cu: the same butterflies, the same order
//of the radix-2 levels, and the same four-step split with the twist after the column transforms (forward) or before them
//(inverse). A row transform is one row here instead of one thread block, and the column transforms are done on blocks of
//NTT_COLUMN_BLOCK columns, copied to a contiguous buffer, so that the butterflies of all columns of a block share a twiddle.

#include "ntt_core.h"
#include "arith32.h"
#include <string.h>

static uint64_t pow_mod(uint64_t a, uint64_t e, uint64_t p)
{
	uint64_t r = 1;
	a %= p;
	while (e != 0)
	{
		if (e & 1)
		{
			r = r * a % p;
		}
		a = a * a % p;
		e >>= 1;
	}
	return r;
}

static uint32_t to_montgomery(uint64_t x, uint32_t p)
{
	return (uint32_t)((x << 32) % p);
}

void Ntt_prime::Init(uint32_t prime)
{
	p = prime;
	uint32_t x = 1;
	for (int i = 0; i < 5; i++)
	{
		x *= 2 - p * x; //Newton iteration for 1 / p mod 2^32
	}
	p_inv = 0u - x;
	generator = 3;
	while (pow_mod(generator, (p - 1) / 2, p) != p - 1)
	{
		generator++;
	}
	two32 = (uint32_t)((1ULL << 32) % p);
	two32_shoup = shoup_factor(two32, p);

	const uint64_t T_max = 1ULL << NTT_MAX_LOG_LENGTH;
	const uint64_t omega = pow_mod(generator, (p - 1) / T_max, p);
	const uint64_t omega_inv = pow_mod(omega, p - 2, p);

	const size_t row = (size_t)1 << NTT_ROW_LOG_LENGTH;
	w.assign(row, 0);
	w_shoup.assign(row, 0);
	iw.assign(row, 0);
	iw_shoup.assign(row, 0);
	for (size_t s = 1; s < row; s *= 2)
	{
		const uint64_t root = pow_mod(omega, T_max / (2 * s), p), root_inv = pow_mod(omega_inv, T_max / (2 * s), p);
		uint64_t a = 1, b = 1;
		for (size_t j = 0; j < s; j++)
		{
			w[s + j] = (uint32_t)a;
			w_shoup[s + j] = shoup_factor((uint32_t)a, p);
			iw[s + j] = (uint32_t)b;
			iw_shoup[s + j] = shoup_factor((uint32_t)b, p);
			a = a * root % p;
			b = b * root_inv % p;
		}
	}

	twist_low.resize(2048);
	itwist_low.resize(2048);
	uint64_t a = 1, b = 1;
	for (size_t t = 0; t < 2048; t++)
	{
		twist_low[t] = to_montgomery(a, p);
		itwist_low[t] = to_montgomery(b, p);
		a = a * omega % p;
		b = b * omega_inv % p;
	}
	const uint64_t step = pow_mod(omega, 2048, p), step_inv = pow_mod(omega_inv, 2048, p);
	twist_high.resize(T_max / 2048);
	itwist_high.resize(T_max / 2048);
	a = 1;
	b = 1;
	for (size_t h = 0; h < T_max / 2048; h++)
	{
		twist_high[h] = to_montgomery(a, p);
		itwist_high[h] = to_montgomery(b, p);
		a = a * step % p;
		b = b * step_inv % p;
	}
}

size_t ntt_scratch_size(int logT)
{
	if (logT <= NTT_ROW_LOG_LENGTH)
	{
		return 0;
	}
	int logR = (logT - NTT_ROW_LOG_LENGTH > 3) ? logT - NTT_ROW_LOG_LENGTH : 3;
	return ((size_t)1 << logR) * NTT_COLUMN_BLOCK;
}

static size_t bit_reverse(size_t i, int bits)
{
	size_t r = 0;
	for (int k = 0; k < bits; k++)
	{
		r = (r << 1) | ((i >> k) & 1);
	}
	return r;
}

//Decimation-in-frequency transform of M = 2^logM contiguous residues, the top level (half size M / 2) first.
static void row_forward(uint32_t* x, int logM, const Ntt_prime& t)
{
	const uint32_t p = t.p;
	const size_t M = (size_t)1 << logM;
	for (size_t s = M / 2; s >= 1; s /= 2)
	{
		const uint32_t* w = t.w.data() + s;
		const uint32_t* wp = t.w_shoup.data() + s;
		for (size_t block = 0; block < M; block += 2 * s)
		{
			uint32_t* u = x + block;
			uint32_t* v = u + s;
			for (size_t j = 0; j < s; j++)
			{
				uint32_t a = u[j], b = v[j];
				u[j] = add_mod(a, b, p);
				v[j] = mul_shoup(sub_mod(a, b, p), w[j], wp[j], p);
			}
		}
	}
}

//Decimation-in-time inverse of row_forward, times M, the bottom level (half size 1) first.
static void row_inverse(uint32_t* x, int logM, const Ntt_prime& t)
{
	const uint32_t p = t.p;
	const size_t M = (size_t)1 << logM;
	for (size_t s = 1; s < M; s *= 2)
	{
		const uint32_t* w = t.iw.data() + s;
		const uint32_t* wp = t.iw_shoup.data() + s;
		for (size_t block = 0; block < M; block += 2 * s)
		{
			uint32_t* u = x + block;
			uint32_t* v = u + s;
			for (size_t j = 0; j < s; j++)
			{
				uint32_t a = u[j], b = mul_shoup(v[j], w[j], wp[j], p);
				u[j] = add_mod(a, b, p);
				v[j] = sub_mod(a, b, p);
			}
		}
	}
}

//The twist factor omega_T^E, with E scaled to Tmax, in Montgomery form.
static inline uint32_t twist_factor(const uint32_t* high, const uint32_t* low, size_t E, uint32_t p, uint32_t p_inv)
{
	return mul_mont(high[E >> 11], low[E & 2047], p, p_inv);
}

//Multiply the rows of a block of columns by the (inverse) twist factors: row i, column l by omega_T^(l bitreverse(i)).
static void twist_block(uint32_t* buf, size_t l0, int logT, int logR, const uint32_t* high, const uint32_t* low, const Ntt_prime& t)
{
	const uint32_t p = t.p, p_inv = t.p_inv;
	const int shift = NTT_MAX_LOG_LENGTH - logT;
	const size_t R = (size_t)1 << logR;
	for (size_t i = 0; i < R; i++)
	{
		const size_t k1 = bit_reverse(i, logR);
		uint32_t* row = buf + i * NTT_COLUMN_BLOCK;
		for (size_t c = 0; c < NTT_COLUMN_BLOCK; c++)
		{
			size_t E = ((l0 + c) * k1) << shift;
			row[c] = mul_mont(row[c], twist_factor(high, low, E, p, p_inv), p, p_inv);
		}
	}
}

//The forward column transforms of a four-step transform, and the twist.
static void columns_forward(uint32_t* x, int logT, int logR, const Ntt_prime& t, uint32_t* buf)
{
	const uint32_t p = t.p;
	const size_t R = (size_t)1 << logR, C = (size_t)1 << (logT - logR);
	const size_t B = NTT_COLUMN_BLOCK;
	for (size_t l0 = 0; l0 < C; l0 += B)
	{
		for (size_t i = 0; i < R; i++)
		{
			memcpy(buf + i * B, x + i * C + l0, B * sizeof(uint32_t));
		}
		for (size_t s = R / 2; s >= 1; s /= 2)
		{
			for (size_t block = 0; block < R; block += 2 * s)
			{
				for (size_t j = 0; j < s; j++)
				{
					const uint32_t w = t.w[s + j], wp = t.w_shoup[s + j];
					uint32_t* u = buf + (block + j) * B;
					uint32_t* v = u + s * B;
					for (size_t c = 0; c < B; c++)
					{
						uint32_t a = u[c], b = v[c];
						u[c] = add_mod(a, b, p);
						v[c] = mul_shoup(sub_mod(a, b, p), w, wp, p);
					}
				}
			}
		}
		twist_block(buf, l0, logT, logR, t.twist_high.data(), t.twist_low.data(), t);
		for (size_t i = 0; i < R; i++)
		{
			memcpy(x + i * C + l0, buf + i * B, B * sizeof(uint32_t));
		}
	}
}

//The inverse twist, and the inverse column transforms of a four-step transform.
static void columns_inverse(uint32_t* x, int logT, int logR, const Ntt_prime& t, uint32_t* buf)
{
	const uint32_t p = t.p;
	const size_t R = (size_t)1 << logR, C = (size_t)1 << (logT - logR);
	const size_t B = NTT_COLUMN_BLOCK;
	for (size_t l0 = 0; l0 < C; l0 += B)
	{
		for (size_t i = 0; i < R; i++)
		{
			memcpy(buf + i * B, x + i * C + l0, B * sizeof(uint32_t));
		}
		twist_block(buf, l0, logT, logR, t.itwist_high.data(), t.itwist_low.data(), t);
		for (size_t s = 1; s < R; s *= 2)
		{
			for (size_t block = 0; block < R; block += 2 * s)
			{
				for (size_t j = 0; j < s; j++)
				{
					const uint32_t w = t.iw[s + j], wp = t.iw_shoup[s + j];
					uint32_t* u = buf + (block + j) * B;
					uint32_t* v = u + s * B;
					for (size_t c = 0; c < B; c++)
					{
						uint32_t a = u[c], b = mul_shoup(v[c], w, wp, p);
						u[c] = add_mod(a, b, p);
						v[c] = sub_mod(a, b, p);
					}
				}
			}
		}
		for (size_t i = 0; i < R; i++)
		{
			memcpy(x + i * C + l0, buf + i * B, B * sizeof(uint32_t));
		}
	}
}

void ntt_forward(uint32_t* x, int logT, const Ntt_prime& prime, uint32_t* scratch)
{
	if (logT <= NTT_ROW_LOG_LENGTH)
	{
		row_forward(x, logT, prime);
		return;
	}
	const int logR = (logT - NTT_ROW_LOG_LENGTH > 3) ? logT - NTT_ROW_LOG_LENGTH : 3;
	const size_t R = (size_t)1 << logR, C = (size_t)1 << (logT - logR);
	columns_forward(x, logT, logR, prime, scratch);
	for (size_t i = 0; i < R; i++)
	{
		row_forward(x + i * C, logT - logR, prime);
	}
}

void ntt_inverse(uint32_t* x, int logT, const Ntt_prime& prime, uint32_t* scratch)
{
	if (logT <= NTT_ROW_LOG_LENGTH)
	{
		row_inverse(x, logT, prime);
		return;
	}
	const int logR = (logT - NTT_ROW_LOG_LENGTH > 3) ? logT - NTT_ROW_LOG_LENGTH : 3;
	const size_t R = (size_t)1 << logR, C = (size_t)1 << (logT - logR);
	for (size_t i = 0; i < R; i++)
	{
		row_inverse(x + i * C, logT - logR, prime);
	}
	columns_inverse(x, logT, logR, prime, scratch);
}

void ntt_pointwise(uint32_t* x, const uint32_t* y, size_t n, const Ntt_prime& prime)
{
	const uint32_t p = prime.p, p_inv = prime.p_inv;
	for (size_t i = 0; i < n; i++)
	{
		x[i] = mul_mont(x[i], y[i], p, p_inv);
	}
}
