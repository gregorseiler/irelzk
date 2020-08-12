#include <stdint.h>
#include <immintrin.h>
#include "params.h"
#include "rounding.h"

void power2round_avx(int32_t a1[N], int32_t a0[N], const int32_t a[N])
{
  int i;
  __m256i f,f0,f1;
  const __m256i mask = _mm256_set1_epi32(-(1 << D));
  const __m256i half = _mm256_set1_epi32((1 << (D-1)) - 1);

  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256((__m256i *)&a[8*i]);
    f1 = _mm256_add_epi32(f,half);
    f0 = _mm256_and_si256(f1,mask);
    f1 = _mm256_srai_epi32(f1,D);
    f0 = _mm256_sub_epi32(f,f0);
    _mm256_store_si256((__m256i *)&a1[8*i], f1);
    _mm256_store_si256((__m256i *)&a0[8*i], f0);
  }
}

void decompose_avx(int32_t a1[N], int32_t a0[N], const int32_t a[N])
{
  int i;
  __m256i f,f0,f1,t0,t1;
  const __m256i half = _mm256_set1_epi32(1 << 17);
  const __m256i mask = _mm256_set1_epi32(-(1 << 18));
  const __m256i max = _mm256_set1_epi32(2048);

  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256((__m256i *)&a[8*i]);
    t0 = _mm256_srai_epi32(f,12);
    t1 = _mm256_srai_epi32(f,24);
    t0 = _mm256_add_epi32(t0,f);
    t1 = _mm256_add_epi32(t1,half);
    t0 = _mm256_add_epi32(t0,t1);
    f1 = _mm256_srai_epi32(t0,18);

    t0 = _mm256_and_si256(t0,mask);
    t1 = _mm256_slli_epi32(f1,6);
    f0 = _mm256_sub_epi32(f,t0);
    f0 = _mm256_add_epi32(f0,t1);

    t0 = _mm256_cmpeq_epi32(f1,max);
    f1 = _mm256_xor_si256(f1,t0);
    f1 = _mm256_sub_epi32(f1,t0);
    f0 = _mm256_add_epi32(f0,t0);

    _mm256_store_si256((__m256i *)&a1[8*i],f1);
    _mm256_store_si256((__m256i *)&a0[8*i],f0);
  }
}

void makehint_avx(int32_t h[N], const int32_t a1[N], const int32_t a0[N])
{
  int i;
  __m256i f0, f1, g0, g1, g2;
  const __m256i blo = _mm256_set1_epi32(-GAMMA2);
  const __m256i bhi = _mm256_set1_epi32(GAMMA2);
  const __m256i min = _mm256_set1_epi32(-2048);

  for(i=0;i<N/8;i++) {
    f0 = _mm256_load_si256((__m256i *)&a0[8*i]);
    f1 = _mm256_load_si256((__m256i *)&a1[8*i]);
    g0 = _mm256_cmpgt_epi32(blo,f0);
    g1 = _mm256_cmpgt_epi32(f0,bhi);
    g0 = _mm256_or_si256(g0,g1);
    g1 = _mm256_cmpeq_epi32(blo,f0);
    g2 = _mm256_cmpeq_epi32(f1,min);
    g1 = _mm256_andnot_si256(g2,g1);
    g0 = _mm256_or_si256(g0,g1);
    g0 = _mm256_sign_epi32(g0,g0);
    _mm256_store_si256((__m256i *)&h[8*i],g0);
  }
}

void usehint_avx(int32_t b1[N], const int32_t a[N], const int32_t hint[N]) {
  int i;
  __attribute__((aligned(32)))
  int32_t a0[N];
  __m256i f, g, h;
  const __m256i off = _mm256_set1_epi32(2048);
  const __m256i mask = _mm256_set1_epi32(4095);

  decompose_avx(b1, a0, a);
  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256((__m256i *)&a0[8*i]);
    g = _mm256_load_si256((__m256i *)&b1[8*i]);
    h = _mm256_load_si256((__m256i *)&hint[8*i]);
    h = _mm256_sign_epi32(h,f);
    g = _mm256_add_epi32(g,h);
    g = _mm256_add_epi32(g,off);
    g = _mm256_and_si256(g,mask);
    g = _mm256_sub_epi32(g,off);
    _mm256_store_si256((__m256i *)&b1[8*i],g);
  }
}
