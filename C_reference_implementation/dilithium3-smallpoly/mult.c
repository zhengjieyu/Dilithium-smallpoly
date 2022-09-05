/*
Parallel version
Corresponding Algorithm 9
*/

#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "mult.h"

void prepare_s1_table(uint64_t s_table[4*N], const polyvecl *s1)
{
	uint32_t k,j;
    uint64_t temp , temp2;
	uint64_t mask_s = 0x8040201008;
   

	for(k=0; k<N; k++){
        //s_table[k+N] = 0;

        for(j=0; j<L; j++)
		{
			temp = (uint64_t)(ETA + s1->vec[j].coeffs[k]);
			temp2 = (uint64_t)(ETA - s1->vec[j].coeffs[k]);
			s_table[k+N] = (s_table[k+N]<<9) | (temp);
			s_table[k] = (s_table[k]<<9) | (temp2);
		}
     
        s_table[k+3*N] = mask_s - s_table[k+N];
        s_table[k+2*N] = mask_s - s_table[k];

    }

}

void evaluate_cs1(polyvecl *z, const poly *c, const uint64_t s_table[4*N])
{
    uint32_t i,j;
    uint64_t answer[N]={0};
	uint64_t temp;

    for( i = 0 ; i < N ; i++){
        if(c->coeffs[i] == 1){
            for(j = 0 ; j < N ; j++){
                answer[j] += s_table[j-i+N];
            }
        }
        else if(c->coeffs[i] == -1){
            for(j = 0 ; j < N ; j++){
                answer[j] += s_table[j-i+3*N];
            }
        }
    }
    
    for(i=0; i<N; i++)
	{
        temp = answer[i];

        for(j=0; j<L; j++)
		{
			z->vec[L-1-j].coeffs[i] = ((int32_t)(temp & 0x1FF)-TAU*ETA);//(t mod M) - rU (mod q)
			temp >>= 9;//t = t/M

		}

	}

}
        


void prepare_s2_table(uint64_t s2_table[4*N],const polyveck *s2)
{
	uint32_t k,j;
    uint64_t temp , temp2;
    uint64_t mask_s = 0x1008040201008;


	for(k=0; k<N; k++){
        s2_table[k+N] = 0;

        for(j=0; j<K; j++)
		{
			temp = (uint64_t)(ETA + s2->vec[j].coeffs[k]);
			temp2 = (uint64_t)(ETA - s2->vec[j].coeffs[k]);
			s2_table[k+N] = (s2_table[k+N]<<9) | (temp);
			s2_table[k] = (s2_table[k]<<9) | (temp2);
		}
     
        s2_table[k+3*N] = mask_s - s2_table[k+N];
        s2_table[k+2*N] = mask_s - s2_table[k];

    }

}

void evaluate_cs2(
		polyveck *z, 
		poly *c, const uint64_t s2_table[4*N])
{
    uint32_t i,j;
    uint64_t answer0[N]={0};
	uint64_t temp;

    for( i = 0 ; i < N ; i++){
        if(c->coeffs[i] == 1){
            for(j = 0 ; j < N ; j++){
                answer0[j] += s2_table[j-i+N];

            }
        }
        else if(c->coeffs[i] == -1){
            for(j = 0 ; j < N ; j++){
                answer0[j] += s2_table[j-i+3*N];

            }
        }
    }
    
    for(i=0; i<N; i++)
	{
        temp = answer0[i];

        for(j=0; j<K; j++)
		{
			z->vec[K-1-j].coeffs[i] = ((int32_t)(temp & 0x1FF)-TAU*ETA);//(t mod M) - rU (mod q)
			temp >>= 9;//t = t/M

		}
	}

}


void prepare_t0_table(uint64_t t00_table[4*N], uint64_t t01_table[4*N],const polyveck *t0)
{
	uint32_t k,j;
    uint64_t temp , temp2;
	uint64_t mask_t00, mask_t01;	
	mask_t00 = 0x8000100002000;
	mask_t01 = 0x8000100002000;

	for(k=0; k<N; k++){
        t00_table[k+N] = 0;
        for(j=0; j<3; j++)
		{
			temp = (uint64_t)(0x1000 + t0->vec[j].coeffs[k]);
			temp2 = (uint64_t)(0x1000 - t0->vec[j].coeffs[k]);
			t00_table[k+N] = (t00_table[k+N]<<19) | (temp);
			t00_table[k] = (t00_table[k]<<19) | (temp2);
		}
     
        t00_table[k+3*N] = mask_t00 - t00_table[k+N];
        t00_table[k+2*N] = mask_t00 - t00_table[k];

        t01_table[k+N] = 0;
        for(j=3; j<6; j++)
		{
			temp = (uint64_t)(0x1000 + t0->vec[j].coeffs[k]);
			temp2 = (uint64_t)(0x1000 - t0->vec[j].coeffs[k]);
			t01_table[k+N] = (t01_table[k+N]<<19) | (temp);
			t01_table[k] = (t01_table[k]<<19) | (temp2);
		}
     
        t01_table[k+3*N] = mask_t01 - t01_table[k+N];
        t01_table[k+2*N] = mask_t01 - t01_table[k];

    }

}

void evaluate_ct0(polyveck *z,const poly *c, 
                    const uint64_t t00_table[4*N], const uint64_t t01_table[4*N])
{
    uint32_t i,j;
    uint64_t answer0[N]={0},answer1[N]={0};
	uint64_t temp;

    for( i = 0 ; i < N ; i++){
        if(c->coeffs[i] == 1){
            for(j = 0 ; j < N ; j++){
                answer0[j] += t00_table[j-i+N];
                answer1[j] += t01_table[j-i+N];
            }
        }
        else if(c->coeffs[i] == -1){
            for(j = 0 ; j < N ; j++){
                answer0[j] += t00_table[j-i+3*N];
                answer1[j] += t01_table[j-i+3*N];
            }
        }
    }
    

    for(i=0; i<N; i++)
	{
        temp = answer0[i];
        for(j=0; j<3; j++)
		{
			z->vec[2-j].coeffs[i] = ((int32_t)(temp & 0x7FFFF)-200704);//(t mod M) - rU (mod q)
			temp >>= 19;//t = t/M

		}

        temp = answer1[i];
        for(j=0; j<3; j++)
		{
			z->vec[5-j].coeffs[i] = ((int32_t)(temp & 0x7FFFF)-200704);//(t mod M) - rU (mod q)
			temp >>= 19;//t = t/M
		}
	}
}
        

void prepare_t1_table(uint64_t t10_table[4*N], uint64_t t11_table[4*N],const polyveck *t1)
{
	uint32_t k,j;
    uint64_t temp , temp2;
	uint64_t mask_t10, mask_t11;
	mask_t10 = 0x200010000800;
	mask_t11 = 0x200010000800;

	for(k=0; k<N; k++){
        t10_table[k+N] = 0;
        for(j=0; j<3; j++)
		{
			temp = (uint64_t)(0x0400 + t1->vec[j].coeffs[k]);
			temp2 = (uint64_t)(0x0400 - t1->vec[j].coeffs[k]);
			t10_table[k+N] = (t10_table[k+N]<<17) | (temp);
			t10_table[k] = (t10_table[k]<<17) | (temp2);
		}
     
        t10_table[k+3*N] = mask_t10 - t10_table[k+N];
        t10_table[k+2*N] = mask_t10 - t10_table[k];


        t11_table[k+N] = 0;
        for(j=3; j<6; j++)
		{
			temp = (uint64_t)(0x0400 + t1->vec[j].coeffs[k]);
			temp2 = (uint64_t)(0x0400 - t1->vec[j].coeffs[k]);
			t11_table[k+N] = (t11_table[k+N]<<17) | (temp);
			t11_table[k] = (t11_table[k]<<17) | (temp2);
		}
     
        t11_table[k+3*N] = mask_t11 - t11_table[k+N];
        t11_table[k+2*N] = mask_t11 - t11_table[k];
    }

}

void evaluate_ct1(polyveck *z, const poly *c, 
		const uint64_t t10_table[4*N], const uint64_t t11_table[4*N])
{
    uint32_t i,j;
    uint64_t answer0[N]={0},answer1[N]={0};
	uint64_t temp;

    for( i = 0 ; i < N ; i++){
        if(c->coeffs[i] == 1){
            for(j = 0 ; j < N ; j++){
                answer0[j] += t10_table[j-i+N];
                answer1[j] += t11_table[j-i+N];
            }
        }
        else if(c->coeffs[i] == -1){
            for(j = 0 ; j < N ; j++){
                answer0[j] += t10_table[j-i+3*N];
                answer1[j] += t11_table[j-i+3*N];
            }
        }
    }
    
    for(i=0; i<N; i++)
	{
        temp = answer0[i];
        for(j=0; j<3; j++)
		{
			z->vec[2-j].coeffs[i] = ((int32_t)(temp & 0x1FFFF)-50176);//(t mod M) - rU (mod q)
			temp >>= 17;//t = t/M
		}

        temp = answer1[i]; 
        for(j=0; j<3; j++)
		{
			z->vec[5-j].coeffs[i] = ((int32_t)(temp & 0x1FFFF)-50176);//(t mod M) - rU (mod q)
			temp >>= 17;//t = t/M
		}
	}
}
        