#include <stdint.h>
#include <immintrin.h>
#include "randombytes.h"
#include "aes256ctr.h"
#include "params.h"
#include "consts.h"
#include "ntt.h"
#include "rounding.h"
#include "poly.h"

static unsigned int rej_uniform(int32_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  int32_t t;

  ctr = pos = 0;
  while(ctr < len && pos + 4 <= buflen) {
    t  = buf[pos++];
    t |= (int32_t)buf[pos++] << 8;
    t |= (int32_t)buf[pos++] << 16;
    t |= (int32_t)buf[pos++] << 24;
    t &= (1 << 30)-1;

    if(t < Q)
      r[ctr++] = t;
  }

  return ctr;
}

#define REJ_UNIFORM_BUFLEN ((512 + AES256CTR_BLOCKBYTES - 1)/AES256CTR_BLOCKBYTES*AES256CTR_BLOCKBYTES)
static unsigned int rej_uniform_avx(int32_t * restrict r, const uint8_t buf[REJ_UNIFORM_BUFLEN]) {
  unsigned int ctr, pos;
  uint32_t good;
  __m256i d, tmp;
  const __m256i bound = _mm256_load_si256((__m256i *)&qdata[_8XQ]);
  const __m256i mask  = _mm256_set1_epi32((1 << 30)-1);

  ctr = pos = 0;
  while(pos <= REJ_UNIFORM_BUFLEN-32) {
    d = _mm256_load_si256((__m256i *)&buf[pos]);
    d = _mm256_and_si256(d,mask);
    pos += 32;

    tmp = _mm256_sub_epi32(d,bound);
    good = _mm256_movemask_ps((__m256)tmp);
    if(good == 255) {
      _mm256_storeu_si256((__m256i *)&r[ctr],d);
      ctr += 8;
      continue;
    }

    tmp = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i *)&rejidx[good]));
    d = _mm256_permutevar8x32_epi32(d,tmp);

    _mm256_storeu_si256((__m256i *)&r[ctr],d);
    ctr += _mm_popcnt_u32(good);
  }

  return ctr;
}

#define POLY_UNIFORM_NBLOCKS ((512+AES256CTR_BLOCKBYTES-1)/AES256CTR_BLOCKBYTES)
void poly_uniform_preinit(poly *r, aes256ctr_ctx *state) {
  unsigned int ctr;
  __attribute__((aligned(32)))
  uint8_t buf[POLY_UNIFORM_NBLOCKS*AES256CTR_BLOCKBYTES];

  aes256ctr_squeezeblocks(buf,POLY_UNIFORM_NBLOCKS,state);
  ctr = rej_uniform_avx(r->coeffs,buf);

  while(ctr < N) {
    /* length of buf is divisible by 4; hence, no bytes left */
    aes256ctr_squeezeblocks(buf,1,state);
    ctr += rej_uniform(&r->coeffs[ctr], N - ctr, buf, AES256CTR_BLOCKBYTES);
  }
}

void poly_uniform(poly *r, const uint8_t seed[SYMBYTES], uint16_t nonce) {
  aes256ctr_ctx state;
  aes256ctr_init(&state,seed,nonce);
  poly_uniform_preinit(r,&state);
}

void poly_trinary_preinit(poly *r, aes256ctr_ctx *state) {
  int i;
  __attribute__((aligned(32)))
  uint8_t buf[N/2];
  __m256i f, g, h, mask32, mask4, mask2;
  const __m256i lut = _mm256_set1_epi32(0xA815);

  mask32 = _mm256_cmpeq_epi32(lut,lut);
  mask4 = _mm256_srli_epi32(mask32,28);
  mask2 = _mm256_srli_epi32(mask32,30);

  aes256ctr_squeezeblocks(buf,1,state);

  for(i=0;i<N/16;i++) {
    f = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i *)&buf[8*i]));
    g = _mm256_srli_epi32(f,4);
    f = _mm256_and_si256(f,mask4);
    f = _mm256_srlv_epi32(lut,f);
    g = _mm256_srlv_epi32(lut,g);
    h = _mm256_unpacklo_epi32(f,g);
    g = _mm256_unpackhi_epi32(f,g);
    f = _mm256_permute2x128_si256(h,g,0x20);
    g = _mm256_permute2x128_si256(h,g,0x31);
    f = _mm256_and_si256(f,mask2);
    g = _mm256_and_si256(g,mask2);
    f = _mm256_add_epi32(f,mask32);
    g = _mm256_add_epi32(g,mask32);
    _mm256_store_si256((__m256i *)&r->coeffs[16*i+0],f);
    _mm256_store_si256((__m256i *)&r->coeffs[16*i+8],g);
  }
}

void poly_trinary(poly *r, const uint8_t seed[SYMBYTES], uint16_t nonce) {
  aes256ctr_ctx state;
  aes256ctr_init(&state,seed,nonce);
  poly_trinary_preinit(r,&state);
}

#define POLY_UNIFORM_GAMMA_NBLOCKS ((304+AES256CTR_BLOCKBYTES-1)/AES256CTR_BLOCKBYTES)
void poly_uniform_gamma_preinit(poly *r, aes256ctr_ctx *state) {
  int i, pos = 0;
  __attribute__((aligned(32)))
  uint8_t buf[POLY_UNIFORM_GAMMA_NBLOCKS*AES256CTR_BLOCKBYTES];
  __m256i f;
  const __m256i mask  = _mm256_set1_epi32(0x7FFFF);
  const __m256i min  = _mm256_set1_epi32(-(1 << 18));
  const __m256i idx32 = _mm256_set_epi32(5,2,7,4,1,6,3,0);
  const __m256i idx8  = _mm256_set_epi8(-1,10, 9, 8,-1, 8, 7, 6,
                                         6, 5, 4, 3,-1, 3, 2, 1,
                                        -1, 9, 8, 7, 7, 6, 5, 4,
                                        -1, 4, 3, 2,-1, 2, 1, 0);

  aes256ctr_squeezeblocks(buf,POLY_UNIFORM_GAMMA_NBLOCKS,state);

  for(i=0;i<N/8;i++) {
    f = _mm256_loadu_si256((__m256i *)&buf[pos]);
    pos += 19;
    f = _mm256_permute4x64_epi64(f,0x94);
    f = _mm256_shuffle_epi8(f,idx8);
    f = _mm256_srlv_epi32(f,idx32);
    f = _mm256_and_si256(f,mask);
    f = _mm256_add_epi32(f,min);
    _mm256_store_si256((__m256i *)&r->coeffs[8*i],f);
  }
}

void poly_uniform_gamma(poly *r, const uint8_t seed[SYMBYTES], uint16_t nonce) {
  aes256ctr_ctx state;
  aes256ctr_init(&state,seed,nonce);
  poly_uniform_gamma_preinit(r,&state);
}

void poly_reduce(poly *r) {
  int i;
  __m256i f,t;
  const __m256i mask = _mm256_set1_epi32((1 << 30) - 1);

  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256((__m256i *)&r->coeffs[8*i]);
    t = _mm256_srai_epi32(f,30);
    f = _mm256_and_si256(f,mask);
    f = _mm256_sub_epi32(f,t);
    t = _mm256_slli_epi32(t,18);
    f = _mm256_add_epi32(f,t);
    _mm256_store_si256((__m256i *)&r->coeffs[8*i],f);
  }
}

void poly_csubq(poly *r) {
  int i;
  __m256i f,t;
  const __m256i q = _mm256_load_si256((__m256i *)&qdata[_8XQ]);

  for(i=0;i<N;i+=8) {
    f = _mm256_load_si256((__m256i *)&r->coeffs[i]);
    f = _mm256_sub_epi32(f,q);
    t = _mm256_srai_epi32(f,31);
    t = _mm256_and_si256(t,q);
    f = _mm256_add_epi32(f,t);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f);
  }
}

void poly_caddq(poly *r) {
  int i;
  __m256i f,t;
  const __m256i q = _mm256_load_si256((__m256i *)&qdata[_8XQ]);

  for(i=0;i<N;i+=8) {
    f = _mm256_load_si256((__m256i *)&r->coeffs[i]);
    t = _mm256_srai_epi32(f,31);
    t = _mm256_and_si256(t,q);
    f = _mm256_add_epi32(f,t);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f);
  }
}

void poly_freeze(poly *r) {
  int i;
  __m256i f,t;
  const __m256i q = _mm256_load_si256((__m256i *)&qdata[_8XQ]);
  const __m256i hq = _mm256_srli_epi32(q,1);

  for(i=0;i<N;i+=8) {
    f = _mm256_load_si256((__m256i *)&r->coeffs[i]);
    t = _mm256_cmpgt_epi32(f,hq);
    t = _mm256_and_si256(t,q);
    f = _mm256_sub_epi32(f,t);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f);
  }
}

void poly_neg(poly *r, const poly *a) {
  int i;
  __m256i f;
  const __m256i zero = _mm256_setzero_si256();

  for(i=0;i<N;i+=8) {
    f = _mm256_load_si256((__m256i *)&a->coeffs[i]);
    f = _mm256_sub_epi32(zero,f);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f);
  }
}

void poly_add(poly *r, const poly *a, const poly *b) {
  int i;
  __m256i f,g;

  for(i=0;i<N;i+=8) {
    f = _mm256_load_si256((__m256i *)&a->coeffs[i]);
    g = _mm256_load_si256((__m256i *)&b->coeffs[i]);
    f = _mm256_add_epi32(f,g);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f);
  }

  poly_reduce(r);
}

void poly_sub(poly *r, const poly *a, const poly *b) {
  int i;
  __m256i f,g;

  for(i=0;i<N;i+=8) {
    f = _mm256_load_si256((__m256i *)&a->coeffs[i]);
    g = _mm256_load_si256((__m256i *)&b->coeffs[i]);
    f = _mm256_sub_epi32(f,g);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f);
  }

  poly_reduce(r);
}

void poly_ntt(poly *r) {
  ntt_avx(r->coeffs,qdata);
  poly_reduce(r);
}

void poly_invntt(poly *r) {
  poly_scale_montgomery(r,r,33554432);
  invntt_tomont_avx(r->coeffs,qdata);
}

void poly_invntt_tomont(poly *r) {
  poly_scale_montgomery(r,r,-132153352);
  invntt_tomont_avx(r->coeffs,qdata);
}

void poly_pointwise_montgomery(poly *r, const poly *a, const poly *b)
{
  int i;
  __m256i f0,f1,g0,g1;
  const __m256i q = _mm256_load_si256((__m256i *)&qdata[_8XQ]);
  const __m256i qinv = _mm256_load_si256((__m256i *)&qdata[_8XQINV]);

  for(i=0;i<N;i+=8) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[i]);
    g0 = _mm256_load_si256((__m256i *)&b->coeffs[i]);
    f1 = _mm256_mul_epi32(f0,g0);
    f0 = _mm256_srli_epi64(f0,32);
    g0 = (__m256i)_mm256_movehdup_ps((__m256)g0);
    f0 = _mm256_mul_epi32(f0,g0);
    g1 = _mm256_mul_epi32(f1,qinv);
    g0 = _mm256_mul_epi32(f0,qinv);
    g1 = _mm256_mul_epi32(g1,q);
    g0 = _mm256_mul_epi32(g0,q);
    f1 = _mm256_sub_epi32(f1,g1);
    f0 = _mm256_sub_epi32(f0,g0);
    f1 = (__m256i)_mm256_movehdup_ps((__m256)f1);
    f0 = _mm256_blend_epi32(f1,f0,0xAA);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f0);
  }
}

void poly_scale_montgomery(poly *r, const poly *a, int32_t s)
{
  int i;
  __m256i f0,f1,g0,g1;
  const __m256i q = _mm256_load_si256((__m256i *)&qdata[_8XQ]);
  const __m256i lo = _mm256_set1_epi32(s*QINV);
  const __m256i hi = _mm256_set1_epi32(s);

  for(i=0;i<N;i+=8) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[i]);
    f1 = (__m256i)_mm256_movehdup_ps((__m256)f0);
    g0 = _mm256_mul_epi32(f0,lo);
    g1 = _mm256_mul_epi32(f1,lo);
    f0 = _mm256_mul_epi32(f0,hi);
    f1 = _mm256_mul_epi32(f1,hi);
    g0 = _mm256_mul_epi32(g0,q);
    g1 = _mm256_mul_epi32(g1,q);
    f0 = _mm256_sub_epi32(f0,g0);
    f1 = _mm256_sub_epi32(f1,g1);
    f0 = (__m256i)_mm256_movehdup_ps((__m256)f0);
    f0 = _mm256_blend_epi32(f0,f1,0xAA);
    _mm256_store_si256((__m256i *)&r->coeffs[i],f0);
  }
}

void poly_ntt2(poly r[2]) {
  int i;
  __m256i f0,f1,f2;

  for(i=0;i<N/8;i++) {
    f0 = _mm256_load_si256((__m256i *)&r[0].coeffs[8*i]);
    f1 = _mm256_load_si256((__m256i *)&r[1].coeffs[8*i]);
    f2 = _mm256_sub_epi32(f0,f1);
    f1 = _mm256_add_epi32(f0,f1);
     _mm256_store_si256((__m256i *)&r[0].coeffs[8*i],f1);
     _mm256_store_si256((__m256i *)&r[1].coeffs[8*i],f2);
  }

  poly_pointwise_montgomery(&r[0],&r[0],&twist);
  poly_ntt(&r[0]);
  poly_ntt(&r[1]);
}

void poly_invntt2(poly r[2]) {
  int i;
  __m256i f0,f1,f2;

  invntt_tomont_avx(r[0].coeffs,qdata);
  invntt_tomont_avx(r[1].coeffs,qdata);
  poly_pointwise_montgomery(&r[0],&r[0],&invtwist);

  for(i=0;i<N/8;i++) {
    f0 = _mm256_load_si256((__m256i *)&r[0].coeffs[8*i]);
    f1 = _mm256_load_si256((__m256i *)&r[1].coeffs[8*i]);
    f2 = _mm256_sub_epi32(f0,f1);
    f1 = _mm256_add_epi32(f0,f1);
     _mm256_store_si256((__m256i *)&r[0].coeffs[8*i],f1);
     _mm256_store_si256((__m256i *)&r[1].coeffs[8*i],f2);
  }

  poly_scale_montgomery(&r[0],&r[0],16777216);
  poly_scale_montgomery(&r[1],&r[1],16777216);
}

int poly_chknorm(const poly *a, uint32_t b) {
  int i;
  __m256i f,g,t;
  const __m256i q = _mm256_load_si256((__m256i *)&qdata[_8XQ]);
  const __m256i hq = _mm256_srli_epi32(q,1);
  const __m256i bound = _mm256_set1_epi32(b-1);

  t = _mm256_setzero_si256();
  for(i=0;i<N/8;i++) {
    f = _mm256_load_si256((__m256i *)&a->coeffs[8*i]);
    g = _mm256_cmpgt_epi32(f,hq);
    g = _mm256_and_si256(g,q);
    f = _mm256_sub_epi32(f,g);
    g = _mm256_srai_epi32(f,31);
    f = _mm256_xor_si256(f,g);
    f = _mm256_sub_epi32(f,g);
    f = _mm256_cmpgt_epi32(f,bound);
    t = _mm256_or_si256(t,f);
  }

  return !_mm256_testz_si256(t,t);
}

void poly_rotate(poly *r, const poly *a, int k) {
  int32_t i,j,x;
  poly t;

  for(i=0;i<N;i++) {
    j = i+k;
    x = a->coeffs[i];
    x ^= (-(j&N) >> 31) & (x ^ -x);
    t.coeffs[j&(N-1)] = x;
  }

  *r = t;
}

void poly_sigma(poly *r, const poly *a, int k) {
  int32_t i,j,x;
  poly t;

  j = 0;
  for(i=0;i<N;i++) {
    x = a->coeffs[i];
    x ^= (-(j&N) >> 31) & (x ^ -x);
    t.coeffs[j&(N-1)] = x;
    j += k;
  }

  *r = t;
}

void poly_sigma65_ntt(poly *r, const poly *a) {
  int i;
  __m256i f0,f1,f2,f3;
  for(i=0;i<N/64;i++) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+ 0]);
    f1 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+ 8]);
    f2 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+16]);
    f3 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+24]);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+ 0],f2);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+ 8],f3);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+16],f1);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+24],f0);

    f0 = _mm256_load_si256((__m256i *)&a->coeffs[N/2+32*i+ 0]);
    f1 = _mm256_load_si256((__m256i *)&a->coeffs[N/2+32*i+ 8]);
    f2 = _mm256_load_si256((__m256i *)&a->coeffs[N/2+32*i+16]);
    f3 = _mm256_load_si256((__m256i *)&a->coeffs[N/2+32*i+24]);
    _mm256_store_si256((__m256i *)&r->coeffs[N/2+32*i+ 0],f3);
    _mm256_store_si256((__m256i *)&r->coeffs[N/2+32*i+ 8],f2);
    _mm256_store_si256((__m256i *)&r->coeffs[N/2+32*i+16],f0);
    _mm256_store_si256((__m256i *)&r->coeffs[N/2+32*i+24],f1);
  }
}

void poly_sigma129_ntt(poly *r, const poly *a) {
  int i;
  __m256i f0,f1,f2,f3;
  for(i=0;i<N/32;i++) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+ 0]);
    f1 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+ 8]);
    f2 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+16]);
    f3 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+24]);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+ 0],f1);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+ 8],f0);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+16],f3);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+24],f2);
  }
}

void poly_sigma193_ntt(poly *r, const poly *a) {
  int i;
  __m256i f0,f1,f2,f3;
  for(i=0;i<N/64;i++) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+ 0]);
    f1 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+ 8]);
    f2 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+16]);
    f3 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+24]);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+ 0],f3);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+ 8],f2);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+16],f0);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+24],f1);

    f0 = _mm256_load_si256((__m256i *)&a->coeffs[N/2+32*i+ 0]);
    f1 = _mm256_load_si256((__m256i *)&a->coeffs[N/2+32*i+ 8]);
    f2 = _mm256_load_si256((__m256i *)&a->coeffs[N/2+32*i+16]);
    f3 = _mm256_load_si256((__m256i *)&a->coeffs[N/2+32*i+24]);
    _mm256_store_si256((__m256i *)&r->coeffs[N/2+32*i+ 0],f2);
    _mm256_store_si256((__m256i *)&r->coeffs[N/2+32*i+ 8],f3);
    _mm256_store_si256((__m256i *)&r->coeffs[N/2+32*i+16],f1);
    _mm256_store_si256((__m256i *)&r->coeffs[N/2+32*i+24],f0);
  }
}

void poly_trace65_ntt(poly *r, const poly *a) {
  int i;
  __m256i f0,f1,f2,f3;
  const __m256i mask = _mm256_set1_epi32((1 << 30)  - 1);

  for(i=0;i<N/32;i++) {
    f0 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+ 0]);
    f1 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+ 8]);
    f2 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+16]);
    f3 = _mm256_load_si256((__m256i *)&a->coeffs[32*i+24]);

    f0 = _mm256_add_epi32(f0,f1);
    f2 = _mm256_add_epi32(f2,f3);
    f1 = _mm256_srai_epi32(f0,30);
    f3 = _mm256_srai_epi32(f2,30);
    f0 = _mm256_and_si256(f0,mask);
    f2 = _mm256_and_si256(f2,mask);
    f0 = _mm256_sub_epi32(f0,f1);
    f2 = _mm256_sub_epi32(f2,f3);
    f1 = _mm256_slli_epi32(f1,18);
    f3 = _mm256_slli_epi32(f3,18);
    f0 = _mm256_add_epi32(f0,f1);
    f2 = _mm256_add_epi32(f2,f3);

    f0 = _mm256_add_epi32(f0,f2);
    f1 = _mm256_srai_epi32(f0,30);
    f0 = _mm256_and_si256(f0,mask);
    f0 = _mm256_sub_epi32(f0,f1);
    f1 = _mm256_slli_epi32(f1,18);
    f0 = _mm256_add_epi32(f0,f1);

    _mm256_store_si256((__m256i *)&r->coeffs[32*i+ 0],f0);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+ 8],f0);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+16],f0);
    _mm256_store_si256((__m256i *)&r->coeffs[32*i+24],f0);
  }
}

void poly_power2round(poly *a1, poly *a0, poly *a) {
  poly_reduce(a);
  poly_freeze(a);
  power2round_avx(a1->coeffs,a0->coeffs,a->coeffs);
}

void poly_decompose(poly *a1, poly *a0, poly *a) {
  poly_reduce(a);
  poly_freeze(a);
  decompose_avx(a1->coeffs,a0->coeffs,a->coeffs);
}

void poly_makehint(poly *h, const poly *a1, poly *a0) {
  poly_freeze(a0);
  makehint_avx(h->coeffs,a1->coeffs,a0->coeffs);
}

void poly_usehint(poly *b1, poly *a, const poly *h) {
  poly_reduce(a);
  poly_freeze(a);
  usehint_avx(b1->coeffs,a->coeffs,h->coeffs);
}
