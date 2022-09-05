#ifndef SIGN_NTT_C_H
#define SIGN_NTT_C_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

int ntt_c_crypto_sign_signature(uint8_t *sig, size_t *siglen,
                                const uint8_t *m, size_t mlen,
                                const uint8_t *sk);

int ntt_c_crypto_sign(uint8_t *sm, size_t *smlen,
                      const uint8_t *m, size_t mlen,
                      const uint8_t *sk);

int ntt_c_crypto_sign_verify(const uint8_t *sig, size_t siglen,
                             const uint8_t *m, size_t mlen,
                             const uint8_t *pk);

int ntt_c_crypto_sign_open(uint8_t *m, size_t *mlen,
                           const uint8_t *sm, size_t smlen,
                           const uint8_t *pk);

#endif
