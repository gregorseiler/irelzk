#include <stdint.h>
#include <stdio.h>
#include <immintrin.h>
#include "cpucycles.h"
#include "speed_print.h"
#include "randombytes.h"
#include "main.h"
#include "comm.h"
#include "consts.h"

#define NTESTS 500

int main(void) {
  int i, j, k, ret;
  uint64_t tsc[NTESTS];
  uint8_t rho[SYMBYTES];
  poly A[8][8], u[8], s[16];
  comm t;
  commrnd r;
  proof p;
  poly tmp, s1hat[8];
  aes256ctr_ctx state;
  __m256i f,g;
  const __m256i idx32 = _mm256_set_epi32(0,1,2,3,4,5,6,7);

  randombytes(rho,SYMBYTES);
  aes256ctr_init(&state,rho,0);
  for(i=0;i<16;i++) {
    aes256ctr_select(&state,i);
    poly_trinary_preinit(&s[i], &state);
  }
  for(i=0;i<8;i++) {
    s1hat[i] = s[i];
    poly_ntt(&s1hat[i]);
  }

  randombytes(rho,SYMBYTES);
  aes256ctr_init(&state,rho,0);
  for(i=0;i<8;i++) {
    for(j=0;j<8;j++) {
      aes256ctr_select(&state,(i << 16) + j);
      poly_uniform_preinit(&A[i][j],&state);
      poly_pointwise_montgomery(&tmp,&A[i][j],&s1hat[j]);
      if(j==0) u[i] = tmp;
      else poly_add(&u[i],&u[i],&tmp);
    }
    poly_invntt_tomont(&u[i]);
    poly_add(&u[i],&u[i],&s[8+i]);
  }

  for(i=0;i<8;i++) {
    for(j=i;j<8;j++) {
      for(k=0; k < ((j==i) ? N/16 : N/8); k++) {
        f = _mm256_load_si256((__m256i *)&A[i][j].coeffs[8*k]);
        g = _mm256_load_si256((__m256i *)&A[j][i].coeffs[N-8*(k+1)]);
        f = _mm256_permutevar8x32_epi32(f,idx32);
        g = _mm256_permutevar8x32_epi32(g,idx32);
        _mm256_store_si256((__m256i *)&A[j][i].coeffs[N-8*(k+1)],f);
        _mm256_store_si256((__m256i *)&A[i][j].coeffs[8*k],g);
      }
    }
  }

  main_proof(&p,&t,&r,rho,A,u,s);
  ret = main_proof_verify(&p,&t,rho,A,u);
  if(ret) printf("FAILURE!\n");

  for(i=0;i<NTESTS;i++) {
    tsc[i] = cpucycles();
    main_proof(&p,&t,&r,rho,A,u,s);
  }
  print_results("Mainproto Prover: ",tsc,NTESTS);

  for(i=0;i<NTESTS;i++) {
    tsc[i] = cpucycles();
    main_proof_verify(&p,&t,rho,A,u);
  }
  print_results("Mainproto Verifier: ",tsc,NTESTS);

  return 0;
}
