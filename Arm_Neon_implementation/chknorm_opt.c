#include "chknorm_opt.h"



int chknorml(polyvecl* z, int32_t B, polyvecl* y, poly* cp, polyvecl* s1){
  unsigned int i, j;
  int32_t t;

  for(i = 0; i < L; ++i){
    // poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
    // DBENCH_START();

    if(B > (Q-1)/8)
      return 1;

    for(j = 0; j < N; ++j){
      z->vec[i].coeffs[j] = montgomery_reduce((int64_t)cp->coeffs[j] * s1->vec[i].coeffs[j]);
    }

    invntt_tomont(z->vec[i].coeffs);

    for(j = 0; j < N; ++j){
      z->vec[i].coeffs[j] = z->vec[i].coeffs[j] + y->vec[i].coeffs[j];
      z->vec[i].coeffs[j] = reduce32(z->vec[i].coeffs[j]);

      // chknorm
      t = z->vec[i].coeffs[j] >> 31;
      t = z->vec[i].coeffs[j] - (t & 2*z->vec[i].coeffs[j]);

      if(t >= B) {
        // DBENCH_STOP(*tsample);
        return 1;
      }
    }
    // DBENCH_STOP(*tmul);
  }

  // DBENCH_STOP(*tsample);
  return 0;
}

int chknormk(polyveck* w0, int32_t B, polyveck* h, poly* cp, polyveck* s2){
  unsigned int i, j;
  int32_t t;

  for(i = 0; i < K; ++i){
    // poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
    // DBENCH_START();

    if(B > (Q-1)/8)
      return 1;

    for(j = 0; j < N; ++j){
      h->vec[i].coeffs[j] = montgomery_reduce((int64_t)cp->coeffs[j] * s2->vec[i].coeffs[j]);
    }
    
    invntt_tomont(h->vec[i].coeffs);

    for(j = 0; j < N; ++j){
      w0->vec[i].coeffs[j] = w0->vec[i].coeffs[j] - h->vec[i].coeffs[j];
      w0->vec[i].coeffs[j] = reduce32(w0->vec[i].coeffs[j]);

      // chknorm
      t = w0->vec[i].coeffs[j] >> 31;
      t = w0->vec[i].coeffs[j] - (t & 2*w0->vec[i].coeffs[j]);

      if(t >= B) {
        // DBENCH_STOP(*tsample);
        return 1;
      }
    }
    // DBENCH_STOP(*tmul);
  }

  // DBENCH_STOP(*tsample);
  return 0;
}

int chknormk_hint(polyveck* h, int32_t B, poly* cp, polyveck* t0){
  unsigned int i, j;
  int32_t t;

  for(i = 0; i < K; ++i){
    // poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
    // DBENCH_START();

    if(B > (Q-1)/8)
      return 1;

    for(j = 0; j < N; ++j){
      h->vec[i].coeffs[j] = montgomery_reduce((int64_t)cp->coeffs[j] * t0->vec[i].coeffs[j]);
    }
    
    invntt_tomont(h->vec[i].coeffs);

    for(j = 0; j < N; ++j){
      h->vec[i].coeffs[j] = reduce32(h->vec[i].coeffs[j]);

      // chknorm
      t = h->vec[i].coeffs[j] >> 31;
      t = h->vec[i].coeffs[j] - (t & 2*h->vec[i].coeffs[j]);

      if(t >= B) {
        // DBENCH_STOP(*tsample);
        return 1;
      }
    }
    // DBENCH_STOP(*tmul);
  }

  // DBENCH_STOP(*tsample);
  return 0;
}