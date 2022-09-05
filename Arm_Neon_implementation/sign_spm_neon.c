#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec_ntt_c.h"
#include "poly_ntt_c.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"
#include "mult.h"

/*************************************************
 * Name:        crypto_sign_signature
 *
 * Description: Computes signature.
 *
 * Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
 *              - size_t *siglen: pointer to output length of signature
 *              - uint8_t *m:     pointer to message to be signed
 *              - size_t mlen:    length of message
 *              - uint8_t *sk:    pointer to bit-packed secret key
 *
 * Returns 0 (success)
 **************************************************/
int spm_neon_crypto_sign_signature(uint8_t *sig,
                                   size_t *siglen,
                                   const uint8_t *m,
                                   size_t mlen,
                                   const uint8_t *sk)
{
  unsigned int n;
  uint8_t seedbuf[2 * SEEDBYTES + 3 * CRHBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z;
  polyveck t0, s2, w1, w0, h;
  poly cp;
  keccak_state state;

  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + CRHBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;
  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

  /* Compute CRH(tr, msg) */
  shake256_init(&state);
  shake256_absorb(&state, tr, CRHBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
  crh(rhoprime, key, SEEDBYTES + CRHBYTES);
#endif

#if DILITHIUM_MODE == 2
  uint64_t se_table[3 * N] = {0}; //  for s1[..], s2[..]
  uint64_t answer_se[N] = {0};
  uint64_t t00_table[3 * N] = {0}; // for t0[0], t0[1]
  uint64_t t01_table[3 * N] = {0}; // for t0[2], t0[3]
  uint16_t table_offset[N] = {0};

#elif DILITHIUM_MODE == 3
  uint64_t s_table[3 * N] = {0};
  uint64_t e_table[3 * N] = {0};
  uint64_t t00_table[3 * N] = {0}; // for t0[0], t0[1]
  uint64_t t01_table[3 * N] = {0}; // for t0[2], t0[3]
  uint16_t table_offset[N] = {0};

#elif DILITHIUM_MODE == 5
  uint64_t s_table[3 * N] = {0};
  uint64_t e0_table[3 * N] = {0};  // for e[0], e[1], e[2], e[3]
  uint64_t e1_table[3 * N] = {0};  // for e[4], e[5], e[6], e[7]
  uint64_t t00_table[3 * N] = {0}; // for t0[0], t0[1], t0[2]
  uint64_t t01_table[3 * N] = {0}; // for t0[3], t0[4], t0[5]
  uint64_t t02_table[3 * N] = {0}; // for t0[6], t0[7]
  uint16_t table_offset[N] = {0};

#endif

// prepare s_table, e_table, t0_table;
#if DILITHIUM_MODE == 2
  prepare_se_table(se_table, s1.vec, s2.vec);
  prepare_t0_table(t00_table, t01_table, t0.vec);

#elif DILITHIUM_MODE == 3
  prepare_s_table(s_table, s1.vec);
  prepare_e_table(e_table, s2.vec);
  prepare_t0_table(t00_table, t01_table, t0.vec);

#elif DILITHIUM_MODE == 5
  prepare_s_table(s_table, s1.vec);
  prepare_e_table(e0_table, e1_table, s2.vec);
  prepare_t0_table(t00_table, t01_table, t02_table, t0.vec);

#endif

  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho);


rej:
  /* Sample intermediate vector y */
  polyvecl_uniform_gamma1(&y, rhoprime, nonce++);
  z = y;
  polyvecl_ntt_c(&z);

  /* Matrix-vector multiplication */
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont_c(&w1);

  /* Decompose w and call the random oracle */
  polyveck_caddq(&w1);
  polyveck_decompose(&w1, &w0, &w1);
  polyveck_pack_w1(sig, &w1);

  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, sig, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(sig, SEEDBYTES, &state);
  poly_challenge(&cp, sig);


  /* Compute z, reject if it reveals secret */
  // compute c*s1.
  prepare_table_offset(table_offset, &cp);
#if DILITHIUM_MODE == 2
  evaluate_cse(answer_se, &cp, se_table, table_offset);
  recover_cs(z.vec, answer_se);

#else
  evaluate_cs(z.vec, &cp, s_table, table_offset);

#endif


  polyvecl_add(&z, &z, &y);
  polyvecl_reduce(&z);
  if (polyvecl_chknorm(&z, GAMMA1 - BETA))
    goto rej;

/* Check that subtracting cs2 does not change high bits of w and low bits
 * do not reveal secret information */
// compute h:=c*s2.
#if DILITHIUM_MODE == 2
  recover_ce(h.vec, answer_se);

#elif DILITHIUM_MODE == 3
  evaluate_ce(h.vec, &cp, e_table, table_offset);

#elif DILITHIUM_MODE == 5
  evaluate_ce(h.vec, &cp, e0_table, e1_table, table_offset);

#endif

  // polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
  // polyveck_invntt_tomont(&h);

  polyveck_sub(&w0, &w0, &h);
  polyveck_reduce(&w0);
  if (polyveck_chknorm(&w0, GAMMA2 - BETA))
    goto rej;

/* Compute hints for w1 */
// compute c*t0.
#if DILITHIUM_MODE == 2 || DILITHIUM_MODE == 3
  evaluate_ct0(h.vec, &cp, t00_table, t01_table, table_offset);

#elif DILITHIUM_MODE == 5
  evaluate_ct0(h.vec, &cp, t00_table, t01_table, t02_table, table_offset);

#endif


  if (polyveck_chknorm(&h, GAMMA2))
    goto rej;

  polyveck_add(&w0, &w0, &h);
  polyveck_caddq(&w0);
  n = polyveck_make_hint(&h, &w0, &w1);
  if (n > OMEGA)
    goto rej;

  /* Write signature */
  pack_sig(sig, sig, &z, &h);
  *siglen = CRYPTO_BYTES;
  return 0;
}

/*************************************************
 * Name:        crypto_sign
 *
 * Description: Compute signed message.
 *
 * Arguments:   - uint8_t *sm: pointer to output signed message (allocated
 *                             array with CRYPTO_BYTES + mlen bytes),
 *                             can be equal to m
 *              - size_t *smlen: pointer to output length of signed
 *                               message
 *              - const uint8_t *m: pointer to message to be signed
 *              - size_t mlen: length of message
 *              - const uint8_t *sk: pointer to bit-packed secret key
 *
 * Returns 0 (success)
 **************************************************/
int spm_neon_crypto_sign(uint8_t *sm,
                         size_t *smlen,
                         const uint8_t *m,
                         size_t mlen,
                         const uint8_t *sk)
{
  size_t i;

  for (i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  spm_neon_crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
  *smlen += mlen;
  return 0;
}

/*************************************************
 * Name:        crypto_sign_verify
 *
 * Description: Verifies signature.
 *
 * Arguments:   - uint8_t *m: pointer to input signature
 *              - size_t siglen: length of signature
 *              - const uint8_t *m: pointer to message
 *              - size_t mlen: length of message
 *              - const uint8_t *pk: pointer to bit-packed public key
 *
 * Returns 0 if signature could be verified correctly and -1 otherwise
 **************************************************/
int spm_neon_crypto_sign_verify(const uint8_t *sig,
                                size_t siglen,
                                const uint8_t *m,
                                size_t mlen,
                                const uint8_t *pk)
{
  unsigned int i;
  uint8_t buf[K * POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t c2[SEEDBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h;
  keccak_state state;

#if DILITHIUM_MODE == 2 || DILITHIUM_MODE == 3
  uint64_t t10_table[3 * N] = {0}; // for t1[0], t1[1]
  uint64_t t11_table[3 * N] = {0}; // for t1[2], t1[3]
  uint16_t table_offset[N] = {0};

#elif DILITHIUM_MODE == 5
  uint64_t t10_table[3 * N] = {0}; // for t1[0], t1[1], t1[2]
  uint64_t t11_table[3 * N] = {0}; // for t1[3], t1[4], t1[5]
  uint64_t t12_table[3 * N] = {0}; // for t1[6], t1[7]
  uint16_t table_offset[N] = {0};

#endif

  if (siglen != CRYPTO_BYTES)
    return -1;

  unpack_pk(rho, &t1, pk);
  if (unpack_sig(c, &z, &h, sig))
    return -1;
  if (polyvecl_chknorm(&z, GAMMA1 - BETA))
    return -1;

// prepare t1_table[].
#if DILITHIUM_MODE == 2 || DILITHIUM_MODE == 3
  prepare_t1_table(t10_table, t11_table, t1.vec);

#elif DILITHIUM_MODE == 5
  prepare_t1_table(t10_table, t11_table, t12_table, t1.vec);

#endif

  /* Compute CRH(CRH(rho, t1), msg) */
  crh(mu, pk, CRYPTO_PUBLICKEYBYTES);
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  poly_challenge(&cp, c);
  polyvec_matrix_expand(mat, rho);

  polyvecl_ntt_c(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont_c(&w1);

  // evaluate c*t1.
  prepare_table_offset(table_offset, &cp);
#if DILITHIUM_MODE == 2 || DILITHIUM_MODE == 3
  evaluate_ct1(t1.vec, &cp, t10_table, t11_table, table_offset);

#elif DILITHIUM_MODE == 5
  evaluate_ct1(t1.vec, &cp, t10_table, t11_table, t12_table, table_offset);

#endif

  polyveck_sub(&w1, &w1, &t1);

  /* Reconstruct w1 */
  polyveck_caddq(&w1);
  polyveck_use_hint(&w1, &w1, &h);
  polyveck_pack_w1(buf, &w1);

  /* Call random oracle and verify challenge */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, buf, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(c2, SEEDBYTES, &state);
  for (i = 0; i < SEEDBYTES; ++i)
    if (c[i] != c2[i])
      return -1;

  return 0;
}

/*************************************************
 * Name:        crypto_sign_open
 *
 * Description: Verify signed message.
 *
 * Arguments:   - uint8_t *m: pointer to output message (allocated
 *                            array with smlen bytes), can be equal to sm
 *              - size_t *mlen: pointer to output length of message
 *              - const uint8_t *sm: pointer to signed message
 *              - size_t smlen: length of signed message
 *              - const uint8_t *pk: pointer to bit-packed public key
 *
 * Returns 0 if signed message could be verified correctly and -1 otherwise
 **************************************************/
int spm_neon_crypto_sign_open(uint8_t *m,
                              size_t *mlen,
                              const uint8_t *sm,
                              size_t smlen,
                              const uint8_t *pk)
{
  size_t i;

  if (smlen < CRYPTO_BYTES)
    goto badsig;

  *mlen = smlen - CRYPTO_BYTES;
  if (spm_neon_crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
    goto badsig;
  else
  {
    /* All good, copy msg, return 0 */
    for (i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_BYTES + i];
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = -1;
  for (i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}