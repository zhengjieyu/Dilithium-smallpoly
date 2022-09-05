#ifndef MULT_C_H
#define MULT_C_H

#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "reduce.h"

void prepare_se_table_c(
		uint64_t se_table[3 * N],
		const poly s[L], const poly e[L]);
void evaluate_cse_c(
		uint64_t answer_se[N],
		const poly *c, const uint64_t se_table[3 * N]);
void recover_cs_c(poly cs[L], uint64_t answer_se[N]);
void recover_ce_c(poly ce[L], const uint64_t answer_se[N]);

void prepare_t0_table_c(
		uint64_t t00_table[3 * N], uint64_t t01_table[3 * N],
		const poly t0[K]);
void prepare_t1_table_c(
		uint64_t t10_table[3 * N], uint64_t t11_table[3 * N],
		const poly t1[K]);

void evaluate_ct0_c(
		poly target[K],
		const poly *c,
		const uint64_t t00_table[3 * N], const uint64_t t01_table[3 * N]);
void evaluate_ct1_c(
		poly target[K],
		const poly *c,
		const uint64_t t10_table[3 * N], const uint64_t t11_table[3 * N]);

#endif
