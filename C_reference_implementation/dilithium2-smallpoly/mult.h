#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"


void prepare_s1_table(uint64_t s_table[4*N], polyvecl *s1_exp);
void evaluate_cs1(polyvecl *z, const poly *c, const uint64_t s_table[4*N]);
void prepare_s2_table(uint64_t s2_table[4*N], polyveck *s2);
void evaluate_cs2(polyveck *z, const poly *c, const uint64_t s2_table[4*N]);
void prepare_t0_table(uint64_t t00_table[4*N], uint64_t t01_table[4*N], polyveck *t0);
void evaluate_ct0(polyveck *z,const poly *c, 
                    const uint64_t t00_table[4*N], const uint64_t t01_table[4*N]);
void prepare_t1_table(uint64_t t10_table[4*N], uint64_t t11_table[4*N], polyveck *t1);
void evaluate_ct1(polyveck *z, const poly *c, 
		const uint64_t t10_table[4*N], const uint64_t t11_table[4*N]);
