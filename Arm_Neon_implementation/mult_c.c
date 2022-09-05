#include "mult.h"

void prepare_se_table_c(
		uint64_t se_table[3 * N],
		const poly s[L], const poly e[L])
{
	uint32_t j, k;
	uint64_t temp;
	uint64_t mask_se = 0x404040404040404;

	for (k = 0; k < N; k++)
	{
		se_table[k + N] = 0;
		for (j = 0; j < L; j++)
		{
			temp = (uint64_t)(e[j].coeffs[k] + ETA);
			se_table[k + N] = (se_table[k + N] << 8) | (temp);
		}
		for (j = 0; j < L; j++)
		{
			temp = (uint64_t)(s[j].coeffs[k] + ETA);
			se_table[k + N] = (se_table[k + N] << 8) | (temp);
		}

		se_table[k] = mask_se - se_table[k + N];
		se_table[k + 2 * N] = se_table[k];
	}
}

void prepare_t0_table_c(
		uint64_t t00_table[3 * N], uint64_t t01_table[3 * N],
		const poly t0[K])
{
	uint32_t j, k;
	uint64_t temp;
	uint64_t mask_t0 = 0x100002000;

	for (k = 0; k < N; k++)
	{
		t00_table[k + N] = 0; // for t0[0], t0[1]
		for (j = 0; j < K / 2; j++)
		{
			temp = (uint64_t)(t0[j].coeffs[k] + (0x1000));
			t00_table[k + N] = (t00_table[k + N] << 19) | (temp);
		}
		t00_table[k] = mask_t0 - t00_table[k + N];
		t00_table[k + 2 * N] = t00_table[k];

		t01_table[k + N] = 0; // for t0[2], t0[3]
		for (j = K / 2; j < K; j++)
		{
			temp = (uint64_t)(t0[j].coeffs[k] + (0x1000));
			t01_table[k + N] = (t01_table[k + N] << 19) | (temp);
		}
		t01_table[k] = mask_t0 - t01_table[k + N];
		t01_table[k + 2 * N] = t01_table[k];
	}
}

void prepare_t1_table_c(
		uint64_t t10_table[3 * N], uint64_t t11_table[3 * N],
		const poly t1[K])
{
	uint32_t j, k;
	uint64_t temp;
	uint64_t mask_t1 = 0x10000800;

	for (k = 0; k < N; k++)
	{
		t10_table[k + N] = 0; // for t1[0], t1[1]
		for (j = 0; j < K / 2; j++)
		{
			temp = (uint64_t)(t1[j].coeffs[k] + (0x0400));
			t10_table[k + N] = (t10_table[k + N] << 17) | (temp);
		}
		t10_table[k] = mask_t1 - t10_table[k + N];
		t10_table[k + 2 * N] = t10_table[k];

		t11_table[k + N] = 0; // for t1[2], t1[3]
		for (j = K / 2; j < K; j++)
		{
			temp = (uint64_t)(t1[j].coeffs[k] + (0x0400));
			t11_table[k + N] = (t11_table[k + N] << 17) | (temp);
		}
		t11_table[k] = mask_t1 - t11_table[k + N];
		t11_table[k + 2 * N] = t11_table[k];
	}
}

void evaluate_cse_c(
		uint64_t answer_se[N],
		const poly *c, const uint64_t se_table[3 * N])
{
	uint32_t j, k;

	for (k = 0; k < N; k++)
		answer_se[k] = 0;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 1)
			for (j = 0; j < N; j++)
				answer_se[j] += se_table[N + j - k];

		if (c->coeffs[k] == -1)
			for (j = 0; j < N; j++)
				answer_se[j] += se_table[2 * N + j - k];
	}
}

void recover_cs_c(poly cs[L], uint64_t answer_se[N])
{
	uint32_t j, k;
	uint64_t temp;

	for (k = 0; k < N; k++)
	{
		temp = answer_se[k];
		for (j = 0; j < L; j++)
		{
			// cs[L-1-j].coeffs[k] = ((int32_t)(temp & 0xFF)-78);
			cs[L - 1 - j].coeffs[k] = (int32_t)((int8_t)temp - 78); // 
			temp >>= 8;
		}
		answer_se[k] = temp;
	}
}

void recover_ce_c(poly ce[L], const uint64_t answer_se[N])
{
	uint32_t j, k;
	uint32_t temp;

	for (k = 0; k < N; k++)
	{
		temp = answer_se[k];
		for (j = 0; j < L; j++)
		{
			ce[L - 1 - j].coeffs[k] = ((int32_t)(temp & 0xFF) - 78);
			temp >>= 8;
		}
	}
}

void evaluate_ct0_c(
		poly target[K],
		const poly *c,
		const uint64_t t00_table[3 * N], const uint64_t t01_table[3 * N])
{
	uint32_t j, k;

	uint64_t answer0[N] = {0}, answer1[N] = {0};
	uint64_t temp;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 1)
			for (j = 0; j < N; j++)
			{
				answer0[j] += t00_table[N + j - k];
				answer1[j] += t01_table[N + j - k];
			}

		if (c->coeffs[k] == -1)
			for (j = 0; j < N; j++)
			{
				answer0[j] += t00_table[2 * N + j - k];
				answer1[j] += t01_table[2 * N + j - k];
			}
	}

	for (k = 0; k < N; k++)
	{
		temp = answer0[k];
		for (j = 0; j < K / 2; j++)
		{
			target[K / 2 - 1 - j].coeffs[k] = ((int32_t)(temp & 0x7FFFF) - 159744);
			temp >>= 19;
		}

		temp = answer1[k];
		for (j = 0; j < K / 2; j++)
		{
			target[K - 1 - j].coeffs[k] = ((int32_t)(temp & 0x7FFFF) - 159744);
			temp >>= 19;
		}
	}
}

void evaluate_ct1_c(
		poly target[K],
		const poly *c,
		const uint64_t t10_table[3 * N], const uint64_t t11_table[3 * N])
{
	uint32_t j, k;
	int32_t t;

	uint64_t answer0[N] = {0}, answer1[N] = {0};
	uint64_t temp;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 1)
			for (j = 0; j < N; j++)
			{
				answer0[j] += t10_table[N + j - k];
				answer1[j] += t11_table[N + j - k];
			}

		if (c->coeffs[k] == -1)
			for (j = 0; j < N; j++)
			{
				answer0[j] += t10_table[2 * N + j - k];
				answer1[j] += t11_table[2 * N + j - k];
			}
	}

	for (k = 0; k < N; k++)
	{
		temp = answer0[k];
		for (j = 0; j < K / 2; j++)
		{
			t = ((int32_t)(temp & 0x1FFFF) - 39936);
			target[K / 2 - 1 - j].coeffs[k] = reduce32(t << D);
			temp >>= 17;
		}

		temp = answer1[k];
		for (j = 0; j < K / 2; j++)
		{
			t = ((int32_t)(temp & 0x1FFFF) - 39936);
			target[K - 1 - j].coeffs[k] = reduce32(t << D);
			temp >>= 17;
		}
	}
}
