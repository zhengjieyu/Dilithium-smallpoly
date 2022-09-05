#ifndef POLY_NTT_C_H
#define POLY_NTT_C_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

void poly_ntt_c(poly *a);

void poly_invntt_tomont_c(poly *a);

#endif
