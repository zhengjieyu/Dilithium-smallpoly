#include "mult.h"
#include "mult_params.h"
#include <string.h>
#include <arm_neon.h>

extern void cal_table(int32_t *coeffs, int32_t add_num, uint64_t *table, int32_t offset);
extern void st_table(uint64_t mask, uint64_t *table);
extern void arrayAcc_arm(uint64_t *array1, uint64_t *array2);
extern void copy_64x8(uint64_t *array1, uint64_t *array2);
extern void recover_cs_arm(uint64_t *temp, int32_t *coeffs, uint32_t and_num, uint32_t sub_num);
extern void evaluate_ct0_arm(uint64_t *temp, int32_t *coeffs, uint32_t and_num, uint32_t sub_num);
extern void evaluate_ct1_arm(uint64_t *temp, int32_t *coeffs, uint32_t and_num, uint32_t sub_num, uint32_t negative_Q);
void prepare_table_offset(uint16_t table_offset[N], const poly *c)
{
	int k;
	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 1)
		{
			table_offset[k] = N - k;
		}
		else if (c->coeffs[k] == -1)
		{
			table_offset[k] = 2 * N - k;
		}
	}
}

#if DILITHIUM_MODE == 2

void prepare_se_table(
		uint64_t se_table[3 * N],
		const poly s[L], const poly e[L])
{
	uint32_t j, k;
	int32_t add_num = U_s;
	uint64_t mask = gamma_s;
	int32_t offset = M_s;
	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < L; j++)
		{
			cal_table(e[j].coeffs + k, add_num, se_table + k + N, offset);
		}
		for (j = 0; j < L; j++)
		{
			cal_table(s[j].coeffs + k, add_num, se_table + k + N, offset);
		}
		st_table(mask, se_table + k + N);
	}
}

void prepare_t0_table(
		uint64_t t00_table[3 * N], uint64_t t01_table[3 * N],
		const poly t0[K])
{
	uint32_t j, k;
	int32_t add_num = U_t0;
	uint64_t mask = gamma_t0;
	int32_t offset = M_t0;
	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < K / 2; j++)
		{
			cal_table(t0[j].coeffs + k, add_num, t00_table + k + N, offset);
		}
		st_table(mask, t00_table + k + N);

		for (j = K / 2; j < K; j++)
		{
			cal_table(t0[j].coeffs + k, add_num, t01_table + k + N, offset);
		}
		
		st_table(mask, t01_table + k + N);
	}
}

void prepare_t1_table(
		uint64_t t10_table[3 * N], uint64_t t11_table[3 * N],
		const poly t1[K])
{
	uint32_t j, k;
	int32_t add_num = U_t1;
	uint64_t mask = gamma_t1;
	int32_t offset = M_t1;
	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < K / 2; j++)
		{
			cal_table(t1[j].coeffs + k, add_num, t10_table + k + N, offset);
		}
		st_table(mask, t10_table + k + N);

		for (j = K / 2; j < K; j++)
		{
			cal_table(t1[j].coeffs + k, add_num, t11_table + k + N, offset);
		}
		st_table(mask, t11_table + k + N);
	}
}

void evaluate_cse(
		uint64_t answer_se[N],
		const poly *c, const uint64_t se_table[3 * N], uint16_t table_offset[N])
{
	uint32_t j, k;
	memset(answer_se, 0, 8 * N);

	for (k = 0; k < N; k++)
	{
		// Pre-record offsets to reduce if
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer_se, se_table + table_offset[k]);
	}
}

void recover_cs(poly cs[L], uint64_t answer_se[N])
{
	uint32_t j, k;
	uint64_t temp[8];
	uint32_t and_num = and_s;
	uint32_t sub_num = TAU_U_s;

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer_se + k);
		for (j = 0; j < L; j++)
		{
			recover_cs_arm(temp, cs[L - 1 - j].coeffs + k, and_num, sub_num);
		}
		copy_64x8(answer_se + k, temp);
	}
}

void recover_ce(poly ce[L], const uint64_t answer_se[N])
{
	uint32_t j, k;
	uint64_t temp[8];
	uint32_t and_num = and_e;
	uint32_t sub_num = TAU_U_e;

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer_se + k);
		for (j = 0; j < L; j++)
		{
			recover_cs_arm(temp, ce[L - 1 - j].coeffs + k, and_num, sub_num);
		}
	}
}

void evaluate_ct0(
		poly target[K],
		const poly *c,
		const uint64_t t00_table[3 * N], const uint64_t t01_table[3 * N],
		uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer0[N] = {0}, answer1[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_t0;
	uint32_t sub_num = TAU_U_t0;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer0, t00_table + table_offset[k]);
		arrayAcc_arm(answer1, t01_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer0 + k);
		for (j = 0; j < K / 2; j++)
		{
			evaluate_ct0_arm(temp, target[K / 2 - 1 - j].coeffs + k, and_num, sub_num);
		}

		copy_64x8(temp, answer1 + k);
		for (j = 0; j < K / 2; j++)
		{
			evaluate_ct0_arm(temp, target[K - 1 - j].coeffs + k, and_num, sub_num);
		}
	}
}

void evaluate_ct1(
		poly target[K],
		const poly *c,
		const uint64_t t10_table[3 * N], const uint64_t t11_table[3 * N],
		uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer0[N] = {0}, answer1[N] = {0};

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer0, t10_table + table_offset[k]);
		arrayAcc_arm(answer1, t11_table + table_offset[k]);
	}

	uint64_t temp[8];
	uint32_t and_num = and_t1;
	uint32_t sub_num = TAU_U_t1;
	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer0 + k);
		for (j = 0; j < K / 2; j++)
		{
			evaluate_ct1_arm(temp, target[K / 2 - 1 - j].coeffs + k, and_num, sub_num, -Q);
		}

		copy_64x8(temp, answer1 + k);
		for (j = 0; j < K / 2; j++)
		{
			evaluate_ct1_arm(temp, target[K - 1 - j].coeffs + k, and_num, sub_num, -Q);
		}
	}
}

#elif DILITHIUM_MODE == 3

void prepare_s_table(uint64_t s_table[3 * N], const poly s[L])
{
	uint32_t j, k;
	int32_t add_num = U_s;
	uint64_t mask_s = gamma_s;
	int32_t offset = M_s;

	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < L; j++)
		{
			cal_table(s[j].coeffs + k, add_num, s_table + k + N, offset);
		}
		st_table(mask_s, s_table + k + N);
	}
}

void prepare_e_table(uint64_t e_table[3 * N], const poly e[K])
{
	uint32_t j, k;
	int32_t add_num = U_e;
	uint64_t mask_e = gamma_e;
	int32_t offset = M_e;

	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < K; j++)
		{
			cal_table(e[j].coeffs + k, add_num, e_table + k + N, offset);
		}
		st_table(mask_e, e_table + k + N);
	}
}

void prepare_t0_table(
		uint64_t t00_table[3 * N], uint64_t t01_table[3 * N],
		const poly t0[K])
{
	uint32_t j, k;
	int32_t add_num = U_t0;
	uint64_t mask_t0 = gamma_t0;
	int32_t offset = M_t0;

	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < K / 2; j++)
		{
			cal_table(t0[j].coeffs + k, add_num, t00_table + k + N, offset);
		}
		
		st_table(mask_t0, t00_table + k + N);

		
		for (j = K / 2; j < K; j++)
		{
			
			cal_table(t0[j].coeffs + k, add_num, t01_table + k + N, offset);
		}
		
		st_table(mask_t0, t01_table + k + N);
	}
}

void prepare_t1_table(
		uint64_t t10_table[3 * N], uint64_t t11_table[3 * N],
		const poly t1[K])
{
	uint32_t j, k;
	int32_t add_num = U_t1;
	uint64_t mask_t1 = gamma_t1;
	int32_t offset = M_t1;

	for (k = 0; k < N; k += 8)
	{
		t10_table[k + N] = 0; // for t1[0], t1[1], t1[2]
		for (j = 0; j < K / 2; j++)
		{
			cal_table(t1[j].coeffs + k, add_num, t10_table + k + N, offset);
		}
		
		st_table(mask_t1, t10_table + k + N);

		t11_table[k + N] = 0; // for t1[3], t1[4], t1[5]
		for (j = K / 2; j < K; j++)
		{
			cal_table(t1[j].coeffs + k, add_num, t11_table + k + N, offset);
		}
		
		st_table(mask_t1, t11_table + k + N);
	}
}

void evaluate_cs(
		poly target[L],
		const poly *c, const uint64_t s_table[3 * N],
		const uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_s;
	uint32_t sub_num = TAU_U_s;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer, s_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		// temp = answer[k];
		copy_64x8(temp, answer + k);
		for (j = 0; j < L; j++)
		{
			recover_cs_arm(temp, target[L - 1 - j].coeffs + k, and_num, sub_num);
		}
	}
}

void evaluate_ce(
		poly target[K],
		const poly *c, const uint64_t e_table[3 * N],
		const uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_e;
	uint32_t sub_num = TAU_U_e;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer, e_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer + k);
		for (j = 0; j < K; j++)
		{
			recover_cs_arm(temp, target[K - 1 - j].coeffs + k, and_num, sub_num);
		}
	}
}

void evaluate_ct0(
		poly target[K],
		const poly *c,
		const uint64_t t00_table[3 * N], const uint64_t t01_table[3 * N],
		const uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer0[N] = {0}, answer1[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_t0;
	uint32_t sub_num = TAU_U_t0;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer0, t00_table + table_offset[k]);
		arrayAcc_arm(answer1, t01_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer0 + k);
		for (j = 0; j < K / 2; j++)
		{
			evaluate_ct0_arm(temp, target[K / 2 - 1 - j].coeffs + k, and_num, sub_num);
		}

		copy_64x8(temp, answer1 + k);
		for (j = 0; j < K / 2; j++)
		{
			evaluate_ct0_arm(temp, target[K - 1 - j].coeffs + k, and_num, sub_num);
		}
	}
}

extern void evaluate_ct1_arm(uint64_t *temp, int32_t *coeffs, uint32_t and_num, uint32_t sub_num, uint32_t negative_Q);
void evaluate_ct1(
		poly target[K],
		const poly *c,
		const uint64_t t10_table[3 * N], const uint64_t t11_table[3 * N],
		const uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer0[N] = {0}, answer1[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_t1;
	uint32_t sub_num = TAU_U_t1;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer0, t10_table + table_offset[k]);
		arrayAcc_arm(answer1, t11_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer0 + k);
		for (j = 0; j < K / 2; j++)
		{
			evaluate_ct1_arm(temp, target[K / 2 - 1 - j].coeffs + k, and_num, sub_num, -Q);
		}

		copy_64x8(temp, answer1 + k);
		for (j = 0; j < K / 2; j++)
		{
			evaluate_ct1_arm(temp, target[K - 1 - j].coeffs + k, and_num, sub_num, -Q);
		}
	}
}

#elif DILITHIUM_MODE == 5

void prepare_s_table(uint64_t s_table[3 * N], const poly s[L])
{
	uint32_t j, k;
	int32_t add_num = U_s;
	uint64_t mask_s = gamma_s;
	int32_t offset = M_s;

	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < L; j++)
		{
			cal_table(s[j].coeffs + k, add_num, s_table + k + N, offset);
		}
		st_table(mask_s, s_table + k + N);
	}
}

void prepare_e_table(
		uint64_t e0_table[3 * N], uint64_t e1_table[3 * N],
		const poly e[K])
{
	uint32_t j, k;
	int32_t add_num = U_e;
	uint64_t mask_e = gamma_e;
	int32_t offset = M_e;

	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < K / 2; j++)
		{
			cal_table(e[j].coeffs + k, add_num, e0_table + k + N, offset);
		}
		st_table(mask_e, e0_table + k + N);

		for (j = K / 2; j < K; j++)
		{
			cal_table(e[j].coeffs + k, add_num, e1_table + k + N, offset);
		}
		st_table(mask_e, e1_table + k + N);
	}
}

void prepare_t0_table(
		uint64_t t00_table[3 * N], uint64_t t01_table[3 * N], uint64_t t02_table[3 * N],
		const poly t0[K])
{
	uint32_t j, k;
	int32_t add_num = U_t0;
	int32_t offset = M_t0;
	uint64_t mask_t00, mask_t01, mask_t02;
	mask_t00 = gamma_t00;
	mask_t01 = gamma_t01;
	mask_t02 = gamma_t02;

	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < 3; j++)
		{
			cal_table(t0[j].coeffs + k, add_num, t00_table + k + N, offset);
		}
		st_table(mask_t00, t00_table + k + N);

		for (j = 3; j < 6; j++)
		{
			cal_table(t0[j].coeffs + k, add_num, t01_table + k + N, offset);
		}
		st_table(mask_t01, t01_table + k + N);

		for (j = 6; j < 8; j++)
		{
			cal_table(t0[j].coeffs + k, add_num, t02_table + k + N, offset);
		}
		st_table(mask_t02, t02_table + k + N);
	}
}

void prepare_t1_table(
		uint64_t t10_table[3 * N], uint64_t t11_table[3 * N], uint64_t t12_table[3 * N],
		const poly t1[K])
{
	uint32_t j, k;
	int32_t add_num = U_t1;
	int32_t offset = M_t1;
	uint64_t mask_t10, mask_t11, mask_t12;
	mask_t10 = gamma_t10;
	mask_t11 = gamma_t11;
	mask_t12 = gamma_t12;

	for (k = 0; k < N; k += 8)
	{
		for (j = 0; j < 3; j++)
		{
			cal_table(t1[j].coeffs + k, add_num, t10_table + k + N, offset);
		}
		st_table(mask_t10, t10_table + k + N);

		for (j = 3; j < 6; j++)
		{
			cal_table(t1[j].coeffs + k, add_num, t11_table + k + N, offset);
		}
		st_table(mask_t11, t11_table + k + N);

		for (j = 6; j < 8; j++)
		{
			cal_table(t1[j].coeffs + k, add_num, t12_table + k + N, offset);
		}
		st_table(mask_t12, t12_table + k + N);
	}
}

void evaluate_cs(
		poly target[L],
		const poly *c, const uint64_t s_table[3 * N],
		const uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_s;
	uint32_t sub_num = TAU_U_s;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer, s_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer + k);
		for (j = 0; j < L; j++)
		{
			recover_cs_arm(temp, target[L - 1 - j].coeffs + k, and_num, sub_num);
		}
	}
}

void evaluate_ce(
		poly target[K],
		const poly *c,
		const uint64_t e0_table[3 * N], const uint64_t e1_table[3 * N],
		const uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer0[N] = {0}, answer1[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_e;
	uint32_t sub_num = TAU_U_e;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer0, e0_table + table_offset[k]);
		arrayAcc_arm(answer1, e1_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer0 + k);
		for (j = 0; j < K / 2; j++)
		{
			recover_cs_arm(temp, target[K / 2 - 1 - j].coeffs + k, and_num, sub_num);
		}

		copy_64x8(temp, answer1 + k);
		for (j = 0; j < K / 2; j++)
		{
			recover_cs_arm(temp, target[K - 1 - j].coeffs + k, and_num, sub_num);
		}
	}
}

extern void evaluate_ct0_arm(uint64_t *temp, int32_t *coeffs, uint32_t and_num, uint32_t sub_num);
void evaluate_ct0(
		poly target[K],
		const poly *c,
		const uint64_t t00_table[3 * N], const uint64_t t01_table[3 * N], const uint64_t t02_table[3 * N],
		const uint16_t table_offset[N])
{
	uint32_t j, k;
	uint64_t answer0[N] = {0}, answer1[N] = {0}, answer2[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_t0;
	uint32_t sub_num = TAU_U_t0;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer0, t00_table + table_offset[k]);
		arrayAcc_arm(answer1, t01_table + table_offset[k]);
		arrayAcc_arm(answer2, t02_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer0 + k);
		for (j = 0; j < 3; j++)
		{
			evaluate_ct0_arm(temp, target[2 - j].coeffs + k, and_num, sub_num);
		}

		copy_64x8(temp, answer1 + k);
		for (j = 0; j < 3; j++)
		{
			evaluate_ct0_arm(temp, target[5 - j].coeffs + k, and_num, sub_num);
		}

		copy_64x8(temp, answer2 + k);
		for (j = 0; j < 2; j++)
		{
			evaluate_ct0_arm(temp, target[7 - j].coeffs + k, and_num, sub_num);
		}
	}
}

extern void evaluate_ct1_arm(uint64_t *temp, int32_t *coeffs, uint32_t and_num, uint32_t sub_num, uint32_t negative_Q);
void evaluate_ct1(
		poly target[K],
		const poly *c,
		const uint64_t t10_table[3 * N], const uint64_t t11_table[3 * N], const uint64_t t12_table[3 * N],
		const uint16_t table_offset[N])
{
	uint32_t j, k;
	int32_t t;
	uint64_t answer0[N] = {0}, answer1[N] = {0}, answer2[N] = {0};
	uint64_t temp[8];
	uint32_t and_num = and_t1;
	uint32_t sub_num = TAU_U_t1;

	for (k = 0; k < N; k++)
	{
		if (c->coeffs[k] == 0)
			continue;
		arrayAcc_arm(answer0, t10_table + table_offset[k]);
		arrayAcc_arm(answer1, t11_table + table_offset[k]);
		arrayAcc_arm(answer2, t12_table + table_offset[k]);
	}

	for (k = 0; k < N; k += 8)
	{
		copy_64x8(temp, answer0 + k);
		for (j = 0; j < 3; j++)
		{
			evaluate_ct1_arm(temp, target[2 - j].coeffs + k, and_num, sub_num, -Q);
		}

		copy_64x8(temp, answer1 + k);
		for (j = 0; j < 3; j++)
		{
			evaluate_ct1_arm(temp, target[5 - j].coeffs + k, and_num, sub_num, -Q);
		}

		copy_64x8(temp, answer2 + k);
		for (j = 0; j < 2; j++)
		{
			evaluate_ct1_arm(temp, target[7 - j].coeffs + k, and_num, sub_num, -Q);
		}
	}
}
#endif