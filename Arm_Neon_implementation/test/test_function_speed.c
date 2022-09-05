#include <stdint.h>
#include "../sign.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../params.h"
#include "cpucycles.h"
#include "speed_print.h"
#include <time.h>
#include "../mult.h"
#include "../mult_c.h"
#include "../packing.h"

int64_t getCpuCycles(void){
  struct timespec time;

  clock_gettime(CLOCK_REALTIME, &time);
  return (int64_t)(time.tv_sec*1e9 + time.tv_nsec);
}

#define MLEN 59
#define NTESTS 1000
uint64_t t[NTESTS];

int main(){
    unsigned int i;
    size_t smlen;
    uint8_t m[MLEN] = {0};
    uint8_t sm[MLEN + CRYPTO_BYTES];
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t seed[CRHBYTES] = {0};
    polyvecl mat[K];
    poly *a = &mat[0].vec[0];
    poly *b = &mat[0].vec[1];
    poly *c = &mat[0].vec[2];

    uint8_t seedbuf[2*SEEDBYTES + 3*CRHBYTES];
    uint8_t *rho, *tr, *key, *mu, *rhoprime;
    rho = seedbuf;
    tr = rho + SEEDBYTES;
    key = tr + CRHBYTES;
    mu = key + SEEDBYTES;
    rhoprime = mu + CRHBYTES;
    polyveck t0;
    uint64_t se_table[3*N] ={0};
    polyvecl s1;
    polyveck s2;

    crypto_sign_keypair(pk, sk);

    unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

    int64_t t1, t2;
    // prepare_se_table_c
    t1 = getCpuCycles();
    for(int i=0;i<NTESTS;i++){
        prepare_se_table_c(se_table, s1.vec, s2.vec);
    }
    t2 = getCpuCycles();
    printf("prepare_se_table_c: %ld\n", t2 - t1);

    // prepare_se_table
    t1 = getCpuCycles();
    for(int i=0;i<NTESTS;i++){
        prepare_se_table(se_table, s1.vec, s2.vec);
    }
    t2 = getCpuCycles();
    printf("prepare_se_table: %ld\n", t2 - t1);

    
    

    return 0;
}