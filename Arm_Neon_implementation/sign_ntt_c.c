#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec_ntt_c.h"
#include "poly_ntt_c.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"

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
int ntt_c_crypto_sign_signature(uint8_t *sig,
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

  // seedbuf = rho || tr || key || mu || rhoprime
  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + CRHBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;
  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

  /* Compute CRH(tr, msg) */
  // mu = CRH(tr || m)
  shake256_init(&state);
  shake256_absorb(&state, tr, CRHBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else 
  crh(rhoprime, key, SEEDBYTES + CRHBYTES); // rhoprime = CRH(key || mu)
#endif

  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho); // mat = ExpandA(rho)
  polyvecl_ntt_c(&s1);
  polyveck_ntt_c(&s2);
  polyveck_ntt_c(&t0);

rej:
  /* Sample intermediate vector y */
  polyvecl_uniform_gamma1(&y, rhoprime, nonce++); // y = ExpandMask(rhoprime, nonce)
  z = y;
  polyvecl_ntt_c(&z);

  /* Matrix-vector multiplication */
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z); // w = Ay
  polyveck_reduce(&w1);
  polyveck_invntt_tomont_c(&w1);

  /* Decompose w and call the random oracle */
  polyveck_caddq(&w1);             
  polyveck_decompose(&w1, &w0, &w1); // w = w1*2^d + w0
  polyveck_pack_w1(sig, &w1);

  // sig = H(mu || w1), cp = SampleInBall(H(mu || w1))
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, sig, K * POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(sig, SEEDBYTES, &state);
  poly_challenge(&cp, sig);
  poly_ntt_c(&cp);

  /* Compute z, reject if it reveals secret */
  polyvecl_pointwise_poly_montgomery(&z, &cp, &s1); // z = cp*s1
  polyvecl_invntt_tomont(&z);
  polyvecl_add(&z, &z, &y); // z = y + cp*s1
  polyvecl_reduce(&z);
  if (polyvecl_chknorm(&z, GAMMA1 - BETA)) // if ||z||∞ >= GAMMA1 - BETA
    goto rej;

  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  polyveck_pointwise_poly_montgomery(&h, &cp, &s2); // h = cp*s2
  polyveck_invntt_tomont(&h);
  polyveck_sub(&w0, &w0, &h); // w0 = w0 - cp*s2
  polyveck_reduce(&w0);
  if (polyveck_chknorm(&w0, GAMMA2 - BETA)) // if ||w0||∞ >= GAMMA2 - BETA
    goto rej;

  /* Compute hints for w1 */
  polyveck_pointwise_poly_montgomery(&h, &cp, &t0); // h = cp*t0
  polyveck_invntt_tomont(&h);
  polyveck_reduce(&h);
  if (polyveck_chknorm(&h, GAMMA2)) // if ||cp*t0||∞ >= GAMMA2
    goto rej;

  polyveck_add(&w0, &w0, &h); // w0 = w0 - cp*s2 + cp*t0
  polyveck_caddq(&w0);
  n = polyveck_make_hint(&h, &w0, &w1); 
  if (n > OMEGA)                        
    goto rej;

  /* Write signature */
  pack_sig(sig, sig, &z, &h); // sig = (H(mu || w1), z, h)
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
int ntt_c_crypto_sign(uint8_t *sm,
                      size_t *smlen,
                      const uint8_t *m,
                      size_t mlen,
                      const uint8_t *sk)
{
  size_t i;

  for (i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  ntt_c_crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
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
int ntt_c_crypto_sign_verify(const uint8_t *sig,
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

  if (siglen != CRYPTO_BYTES)
    return -1;

  unpack_pk(rho, &t1, pk);
  if (unpack_sig(c, &z, &h, sig))
    return -1;
  if (polyvecl_chknorm(&z, GAMMA1 - BETA))
    return -1;

  /* Compute CRH(CRH(rho, t1), msg) */
  // mu = CRH(CRH(rho || t1) || m)
  crh(mu, pk, CRYPTO_PUBLICKEYBYTES); // CRH(rho || t1)
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  poly_challenge(&cp, c);          // cp = SampleInBall(c)
  polyvec_matrix_expand(mat, rho); // mat = ExpandA(rho)

  polyvecl_ntt_c(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z); // w1 = Az

  poly_ntt_c(&cp);
  polyveck_shiftl(&t1); // t1 = t1*2^d
  polyveck_ntt_c(&t1);
  polyveck_pointwise_poly_montgomery(&t1, &cp, &t1); // t1 = cp*t1*2^d

  polyveck_sub(&w1, &w1, &t1); // w1 = Az - c*t1*2^d
  polyveck_reduce(&w1);
  polyveck_invntt_tomont_c(&w1);

  /* Reconstruct w1 */
  polyveck_caddq(&w1);
  polyveck_use_hint(&w1, &w1, &h); // w1 = UseHint(h, Az - c*t1*2^d)
  polyveck_pack_w1(buf, &w1);

  /* Call random oracle and verify challenge */
  // c2 = H(mu || w1)
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
int ntt_c_crypto_sign_open(uint8_t *m,
                           size_t *mlen,
                           const uint8_t *sm,
                           size_t smlen,
                           const uint8_t *pk)
{
  size_t i;

  if (smlen < CRYPTO_BYTES)
    goto badsig;

  *mlen = smlen - CRYPTO_BYTES;
  if (ntt_c_crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
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