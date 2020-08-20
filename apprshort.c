#include <stdint.h>
#include <string.h>
#include "immintrin.h"
#include "aes256ctr.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "apprshort.h"

__attribute__((aligned(32)))
static uint32_t challenge[N][N/32];

#define APPRSHORT_NBLOCKS (N*N/8/AES256CTR_BLOCKBYTES)
void apprshort_challenge_prehash(const uint8_t chash[SYMBYTES]) {
  aes256ctr_ctx state;
  aes256ctr_init(&state,chash,0);
  aes256ctr_squeezeblocks((uint8_t *)challenge,APPRSHORT_NBLOCKS,&state);
}

void apprshort_challenge_mul(poly *r, const poly *a) {
  int i, j, k;
  const __m256i mask = _mm256_set_epi32(128,64,32,16,8,4,2,1);
  const __m256i mask31 = _mm256_set1_epi32((1 << 30) - 1);
  __m256i f0, f1, g0, g1, h[4];

  for(i=0;i<N/32;i++) {
    for(j=0;j<4;j++)
      h[j] = _mm256_setzero_si256();
    for(k=0;k<N;k++) {
      f0 = _mm256_set1_epi32(a->coeffs[k]);
      g0 = _mm256_set1_epi32(challenge[k][i]);
      for(j=0;j<4;j++) {
        g1 = _mm256_and_si256(g0,mask);
        g1 = _mm256_cmpeq_epi32(g1,mask);
        f1 = _mm256_and_si256(f0,g1);
        f1 = _mm256_add_epi32(h[j],f1);
        g1 = _mm256_srai_epi32(f1,30);
        f1 = _mm256_and_si256(f1,mask31);
        f1 = _mm256_sub_epi32(f1,g1);
        g1 = _mm256_slli_epi32(g1,18);
        h[j] = _mm256_add_epi32(f1,g1);
        g0 = _mm256_srli_epi32(g0,8);
      }
    }
    for(j=0;j<4;j++)
      _mm256_store_si256((__m256i *)&r->coeffs[32*i+8*j],h[j]);
  }
}

void apprshort_first(poly *yf, const uint8_t seed[SYMBYTES], uint16_t nonce) {
  int i;
  __m256i f;
  const __m256i mask = _mm256_set1_epi32((1U << 27)-1);
  const __m256i min = _mm256_set1_epi32(-(1U << 26));
  aes256ctr_ctx state;

  aes256ctr_init(&state,seed,nonce);
  aes256ctr_squeezeblocks((uint8_t *)yf,sizeof(poly)/AES256CTR_BLOCKBYTES,&state);
  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256((__m256i *)&yf->coeffs[8*i]);
    f = _mm256_and_si256(f,mask);
    f = _mm256_add_epi32(f,min);
    _mm256_store_si256((__m256i *)&yf->coeffs[8*i],f);
  }
}

int apprshort_last(poly *l, const polyvecm *msg) {
  int i, j, k, ret;
  __attribute__((aligned(32)))
  uint32_t t[8];
  const __m256i mask = _mm256_set_epi32(128,64,32,16,8,4,2,1);
  __m256i f, g, h, acc;

  memset(l,0,sizeof(poly));
  for(i=0;i<N;i++) {
    acc = _mm256_setzero_si256();
    for(j=0;j<N/32;j++) {
      g = _mm256_set1_epi32(challenge[i][j]);
      for(k=0;k<4;k++) {
        f = _mm256_load_si256((__m256i *)&msg->vec[M-4].coeffs[32*j+8*k]);
        h = _mm256_and_si256(g,mask);
        h = _mm256_cmpeq_epi32(h,mask);
        f = _mm256_and_si256(f,h);
        acc = _mm256_add_epi32(acc,f);
        g = _mm256_srli_epi32(g,8);
      }
    }

    _mm256_store_si256((__m256i *)&t,acc);
    for(j=0;j<8;j++)
      l->coeffs[i] += t[j];
  }

  poly_add(l,l,&msg->vec[M-3]);
  ret = poly_chknorm(l,(1U << 26) - 8192);
  return ret;
}

int apprshort_verify(const poly *l) {
  int ret;
  ret = poly_chknorm(l,(1U << 26) - 8192);
  return ret;
}
