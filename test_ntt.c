#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "randombytes.h"
#include "params.h"
#include "poly.h"
#include "cpucycles.h"

static int32_t pow_simple(int32_t a, unsigned int e) {
  int32_t r;
  if(e == 0) return 1;
  else if(e == 1) return a;

  r = pow_simple(a,e/2);
  r = (int64_t)r*r % Q;
  if(e&1) r = (int64_t)r*a % Q;

  return r;
}

static int idx(int i) {
  int r;
  r  = i/32;
  i %= 32;
  r *= 32;
  r += 8*(i%4);
  r += i/4;
  return r;
}

static int32_t zeta[N];
static int32_t zetapow[N][N];

int main(void) {
  int i,j;
  uint64_t t[20], overhead;
  int64_t out[N];
  uint8_t seed[SYMBYTES];
  poly a, b;

  overhead = cpucycles_overhead();
  randombytes(seed,SYMBYTES);

  for(i=0;i<N;i++)
    a.coeffs[i] = 0;
  a.coeffs[1] = 1;
  poly_ntt(&a);
  for(i=0;i<N;i++) {
    zeta[i] = a.coeffs[idx(i)] % Q;
    assert((pow_simple(zeta[i],N) + 1) % Q == 0);
    for(j=0;j<i;j++)
      assert((zeta[j] - zeta[i]) % Q);
  }

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      zetapow[i][j] = pow_simple(zeta[i],j);

  for(i=0;i<N;i++) {
    for(j=0;j<N;j++)
      a.coeffs[j] = 0;
    a.coeffs[i] = 1;
    poly_ntt(&a);
    for(j=0;j<N;j++) {
      a.coeffs[idx(j)] %= Q;
      assert((a.coeffs[idx(j)] - zetapow[j][i]) % Q == 0);
    }
  }

  poly_uniform(&a,seed,0);
  for(i=0;i<N;i++) {
    out[i] = 0;
    for(j=0;j<N;j++)
      out[i] += (int64_t)a.coeffs[j]*zetapow[i][j] % Q;
    out[i] %= Q;
  }
  poly_ntt(&a);
  for(i=0;i<N;i++)
    assert((a.coeffs[idx(i)] - out[i]) % Q == 0);

  poly_uniform(&a,seed,1);
  b = a;
  poly_ntt(&a);
  poly_invntt_tomont(&a);
  poly_scale_montgomery(&a,&a,1);
  for(i=0;i<N;++i)
    assert((a.coeffs[i] - b.coeffs[i]) % Q == 0);

  for(i=0;i<20;i++) {
    t[i] = cpucycles();
    poly_ntt(&a);
  }
  for(i=0;i<19;i++)
    printf("ntt: %2d: %lu\n", i+1, t[i+1] - t[i] - overhead);
  for(i=0;i<20;i++) {
    t[i] = cpucycles();
    poly_invntt_tomont(&a);
  }
  for(i=0;i<19;i++)
    printf("invntt: %2d: %lu\n", i, t[i+1] - t[i] - overhead);

  return 0;
}
