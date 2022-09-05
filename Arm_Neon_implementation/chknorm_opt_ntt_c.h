#ifndef CHKNORM_NTT_C
#define CHKNORM_NTT_C

#include <stdint.h>
#include "params.h"
#include "polyvec_ntt_c.h"
#include "poly_ntt_c.h"
#include "reduce.h"
#include "ntt.h"

int chknorml_c(polyvecl* z, int32_t B, polyvecl* y, poly* cp, polyvecl* s1);

int chknormk_c(polyveck* w0, int32_t B, polyveck* h, poly* cp, polyveck* s2);

int chknormk_hint_c(polyveck* h, int32_t B, poly* cp, polyveck* t0);

#endif