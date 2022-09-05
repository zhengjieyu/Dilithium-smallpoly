#ifndef CHKNORM
#define CHKNORM

#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"
#include "reduce.h"
#include "ntt.h"

int chknorml(polyvecl* z, int32_t B, polyvecl* y, poly* cp, polyvecl* s1);

int chknormk(polyveck* w0, int32_t B, polyveck* h, poly* cp, polyveck* s2);

int chknormk_hint(polyveck* h, int32_t B, poly* cp, polyveck* t0);

#endif