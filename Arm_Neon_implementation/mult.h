#ifndef MULT_H
#define MULT_H

#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "reduce.h"

#if DILITHIUM_MODE == 2
void prepare_se_table(
		uint64_t se_table[3 * N],
		const poly s[L], const poly e[L]);
void prepare_t0_table(
		uint64_t t00_table[3 * N], uint64_t t01_table[3 * N],
		const poly t0[K]);
void prepare_t1_table(
		uint64_t t10_table[3 * N], uint64_t t11_table[3 * N],
		const poly t1[K]);
void prepare_table_offset(uint16_t table_offset[N], const poly *c);

void evaluate_cse(
		uint64_t answer_se[N],
		const poly *c, const uint64_t se_table[3 * N], uint16_t table_offset[N]);
void recover_cs(poly cs[L], uint64_t answer_se[N]);
void recover_ce(poly ce[L], const uint64_t answer_se[N]);
void evaluate_ct0(
		poly target[K],
		const poly *c,
		const uint64_t t00_table[3 * N], const uint64_t t01_table[3 * N], uint16_t table_offset[N]);
void evaluate_ct1(
		poly target[K],
		const poly *c,
		const uint64_t t10_table[3 * N], const uint64_t t11_table[3 * N], uint16_t table_offset[N]);

#elif DILITHIUM_MODE == 3
void prepare_s_table(uint64_t s_table[3 * N], const poly s[L]);
void prepare_e_table(uint64_t e_table[3 * N], const poly e[K]);
void prepare_t0_table(
		uint64_t t00_table[3 * N], uint64_t t01_table[3 * N],
		const poly t0[K]);
void prepare_t1_table(
		uint64_t t10_table[3 * N], uint64_t t11_table[3 * N],
		const poly t1[K]);

void evaluate_cs(
		poly target[L],
		const poly *c, const uint64_t s_table[3 * N],
		const uint16_t table_offset[N]);
void evaluate_ce(
		poly target[K],
		const poly *c, const uint64_t e_table[3 * N],
		const uint16_t table_offset[N]);
void evaluate_ct0(
		poly target[K],
		const poly *c,
		const uint64_t t00_table[3 * N], const uint64_t t01_table[3 * N],
		const uint16_t table_offset[N]);
void evaluate_ct1(
		poly target[K],
		const poly *c,
		const uint64_t t10_table[3 * N], const uint64_t t11_table[3 * N],
		const uint16_t table_offset[N]);

#elif DILITHIUM_MODE == 5
void prepare_s_table(uint64_t s_table[3 * N], const poly s[L]);
void prepare_e_table(uint64_t e0_table[3 * N], uint64_t e1_table[3 * N], const poly e[K]);
void prepare_t0_table(
		uint64_t t00_table[3 * N], uint64_t t01_table[3 * N], uint64_t t02_table[3 * N],
		const poly t0[K]);
void prepare_t1_table(
		uint64_t t10_table[3 * N], uint64_t t11_table[3 * N], uint64_t t12_table[3 * N],
		const poly t1[K]);

void evaluate_cs(
		poly target[L],
		const poly *c, const uint64_t s_table[3 * N],
		const uint16_t table_offset[N]);
void evaluate_ce(
		poly target[K],
		const poly *c,
		const uint64_t e0_table[3 * N], const uint64_t e1_table[3 * N],
		const uint16_t table_offset[N]);
void evaluate_ct0(
		poly target[K],
		const poly *c,
		const uint64_t t00_table[3 * N], const uint64_t t01_table[3 * N], const uint64_t t02_table[3 * N],
		const uint16_t table_offset[N]);

void evaluate_ct1(
		poly target[K],
		const poly *c,
		const uint64_t t10_table[3 * N], const uint64_t t11_table[3 * N], const uint64_t t12_table[3 * N],
		const uint16_t table_offset[N]);

#endif

#endif
