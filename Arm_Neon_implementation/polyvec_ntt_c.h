#ifndef POLYVEC_NTT_C_H
#define POLYVEC_NTT_C_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"

void polyvecl_ntt_c(polyvecl *v);

void polyvecl_invntt_tomont_c(polyvecl *v);

void polyveck_ntt_c(polyveck *v);

void polyveck_invntt_tomont_c(polyveck *v);

#endif
