#ifndef MUL_PARAMS_H
#define MUL_PARAMS_H

#include "params.h"

#if DILITHIUM_MODE == 2

#define gamma_s 0x404040404040404
#define gamma_e 0x404040404040404
#define gamma_t0 0x100002000
#define gamma_t1 0x10000800

#define U_s ETA
#define U_e ETA
#define U_t0 0x1000
#define U_t1 0x0400

#define M_s 8
#define M_e 8
#define M_t0 19
#define M_t1 17

#define and_s 0xFF
#define and_e 0xFF
#define and_t0 0x7FFFF
#define and_t1 0x1FFFF

#define TAU_U_s TAU * U_s
#define TAU_U_e TAU * U_e
#define TAU_U_t0 TAU * U_t0
#define TAU_U_t1 TAU * U_t1

#elif DILITHIUM_MODE == 3

#define gamma_s 0x8040201008
#define gamma_e 0x1008040201008
#define gamma_t0 0x8000100002000
#define gamma_t1 0x200010000800

#define U_s ETA
#define U_e ETA
#define U_t0 0x1000
#define U_t1 0x0400

#define M_s 9
#define M_e 9
#define M_t0 19
#define M_t1 17

#define and_s 0x1FF
#define and_e 0x1FF
#define and_t0 0x7FFFF
#define and_t1 0x1FFFF

#define TAU_U_s TAU * U_s
#define TAU_U_e TAU * U_e
#define TAU_U_t0 TAU * U_t0
#define TAU_U_t1 TAU * U_t1

#elif DILITHIUM_MODE == 5

#define gamma_s 0X100804020100804
#define gamma_e 0x20100804
#define gamma_t00 0x8000100002000
#define gamma_t01 gamma_t00
#define gamma_t02 0x8000100002000
#define gamma_t10 0x200010000800
#define gamma_t11 gamma_t10
#define gamma_t12 0x10000800

#define U_s ETA
#define U_e ETA
#define U_t0 0x1000
#define U_t1 0x0400

#define M_s 9
#define M_e 9
#define M_t0 19
#define M_t1 17

#define and_s 0x1FF
#define and_e 0x1FF
#define and_t0 0x7FFFF
#define and_t1 0x1FFFF

#define TAU_U_s TAU * U_s
#define TAU_U_e TAU * U_e
#define TAU_U_t0 TAU * U_t0
#define TAU_U_t1 TAU * U_t1

#endif

#endif