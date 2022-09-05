#include <stdint.h>
#include "params.h"
#include "reduce.h"
#include <arm_neon.h>

/*************************************************
* Name:        montgomery_reduce
*
* Description: For finite field element a with -2^{31}Q <= a <= Q*2^31,
*              compute r \equiv a*2^{-32} (mod Q) such that -Q < r < Q.
*
* Arguments:   - int64_t: finite field element a
*
* Returns r.
**************************************************/
int32_t montgomery_reduce(int64_t a) {
  int32_t t;

  t = (int32_t)a*QINV;
  t = (a - (int64_t)t*Q) >> 32;
  return t;
}

int32x2_t montgomery_reducex2(int64x2_t ax2){
  int32x2_t tx2;
  int64x2_t tmpx2;
  int32x2_t tmp_ax2;
  int32x2_t QINVx2 = vdup_n_s32(QINV);
  int32x2_t Qx2 = vdup_n_s32(Q);

  tmp_ax2 = vmovn_s64(ax2);
  tx2 = vmul_s32(tmp_ax2, QINVx2);
  tmpx2 = vmull_s32(tx2, Qx2);
  tmpx2 = vsubq_s64(ax2, tmpx2);
  tmpx2 = vshrq_n_s64(tmpx2, 32);
  tx2 = vmovn_s64(tmpx2);

  return tx2;
}

/*************************************************
* Name:        reduce32
*
* Description: For finite field element a with a <= 2^{31} - 2^{22} - 1,
*              compute r \equiv a (mod Q) such that -6283009 <= r <= 6283007.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t reduce32(int32_t a) {
  int32_t t;

  t = (a + (1 << 22)) >> 23;
  t = a - t*Q;
  return t;
}

int32x2_t reduce32x2(int32x2_t a) {
  int32x2_t t;

  // t = (a + (1 << 22)) >> 23;
  // t = a - t*Q;
  int32x2_t num_2p22x2 = vdup_n_s32(1 << 22);
  int32x2_t num_Qx2 = vdup_n_s32(Q);
  t = vadd_s32(a, num_2p22x2);
  t = vshr_n_s32(t, 23);
  t = vmul_s32(t, num_Qx2);
  t = vsub_s32(a, t);
  return t;
}

/*************************************************
* Name:        caddq
*
* Description: Add Q if input coefficient is negative.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t caddq(int32_t a) {
  a += (a >> 31) & Q;
  return a;
}

/*************************************************
* Name:        freeze
*
* Description: For finite field element a, compute standard
*              representative r = a mod^+ Q.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t freeze(int32_t a) {
  a = reduce32(a);
  a = caddq(a);
  return a;
}
