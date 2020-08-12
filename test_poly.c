#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "randombytes.h"
#include "params.h"
#include "consts.h"
#include "poly.h"

static uint8_t bitrev7(uint8_t a) {
  uint8_t t;

  t  = (a & 1) << 6;
  t |= (a & 2) << 4;
  t |= (a & 4) << 2;
  t |= (a & 8) << 0;
  t |= (a & 16) >> 2;
  t |= (a & 32) >> 4;
  t |= (a & 64) >> 6;
  return t;
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

static int idxinv(int i) {
  int r;
  r  = i/32;
  r *= 32;
  i %= 32;
  r += 4*(i%8);
  r += i/8;
  return r;
}

static int lut(int i, int k) {
  int r;
  r = bitrev7(idxinv(i));
  r = (r*k + k/2) % 128;
  r = idx(bitrev7(r));
  return r;
}

int main(void) {
  int i;
  uint8_t seed[SYMBYTES];
  poly f,g,h;

  randombytes(seed,SYMBYTES);

  poly_uniform(&f,seed,0);
  for(i=0;i<N;i++)
    printf("%d, ",f.coeffs[i]);
  printf("\n\n");

  g = f;
  poly_freeze(&g);
  poly_sub(&g,&g,&f);
  poly_freeze(&g);
  for(i=0;i<N;i++)
    printf("%d, ",g.coeffs[i]);
  printf("\n\n");

  poly_trinary(&f,seed,0);
  for(i=0;i<N;i++)
    printf("%d, ",f.coeffs[i]);
  printf("\n\n");

  poly_trinary(&g,seed,0);
  poly_ntt(&f);
  poly_ntt(&g);
  poly_pointwise_montgomery(&f,&f,&g);
  poly_invntt_tomont(&f);
  poly_freeze(&f);
  for(i=0;i<N;i++)
    printf("%d, ",f.coeffs[i]);
  printf("\n\n");

  int max = 0;
  for(i=0;i<N;i++)
    if(abs(f.coeffs[i]) > max)
      max = abs(f.coeffs[i]);
  printf("max = %d; %d\n",max,poly_chknorm(&f,25));

  poly_uniform_gamma(&f,seed,0);
  max = 0;
  for(i=0;i<N;i++)
    if(abs(f.coeffs[i]) > max)
      max = abs(f.coeffs[i]);
  printf("max = %d; %d\n",max,poly_chknorm(&f,(1<<17)));

  poly_uniform(&f,seed,0);
  poly_sigma(&g,&f,65);
  poly_ntt(&f);
  poly_ntt(&g);
  poly_freeze(&f);
  poly_freeze(&g);
  for(i=0;i<N;i++)
    h.coeffs[i] = f.coeffs[lut(i,65)];
  poly_sub(&f,&g,&h);
  for(i=0;i<N;i++)
    printf("%d, ",f.coeffs[idx(i)]);
  printf("\n\n");

  poly_uniform(&f,seed,0);
  poly_sigma(&h,&f,65);
  poly_add(&g,&f,&h);
  poly_sigma(&h,&f,129);
  poly_add(&g,&g,&h);
  poly_sigma(&h,&f,193);
  poly_add(&g,&g,&h);
  poly_ntt(&f);
  poly_ntt(&g);
  poly_trace65_ntt(&f,&f);
  poly_sub(&f,&f,&g);
  poly_freeze(&f);
  for(i=0;i<N;i++)
    printf("%d, ",f.coeffs[idx(i)]);
  printf("\n\n");

  for(i=0;i<N;i++)
    f.coeffs[i] = 0;
  f.coeffs[1] = 1;
  poly_ntt(&f);
  //poly_scale_montgomery(&f,&f,MONTSQ);
  poly_freeze(&f);
  for(i=0;i<N;i++)
    printf("%d, ",f.coeffs[i]);
  printf("\n\n");

  poly a[2];
  for(i=0;i<N;i++)
    a[0].coeffs[i] = a[1].coeffs[i] = 0;
  a[0].coeffs[1] = 1;
  poly_ntt2(a);
  poly_freeze(&a[0]);
  poly_freeze(&a[1]);
  for(i=0;i<N;i++)
    printf("%d, ",a[0].coeffs[i]);
  for(i=0;i<N;i++)
    printf("%d, ",a[1].coeffs[i]);
  printf("\n\n");

  poly_uniform(&f,seed,0);
  poly_decompose(&g,&h,&f);
  for(i=0;i<N;i++)
    if(abs(g.coeffs[i]) > 2048) printf("ERROR\n");
  for(i=0;i<N;i++)
    g.coeffs[i] *= (Q-1)/(1 << 12);
  poly_add(&g,&g,&h);
  poly_sub(&f,&f,&g);
  poly_freeze(&f);
  for(i=0;i<N;i++)
    printf("%d, ",f.coeffs[i]);
  printf("\n\n");
  printf("%d\n", poly_chknorm(&h,(Q-1)/(1 << 13)));

  poly_uniform(&f,seed,0);
  poly_power2round(&g,&h,&f);
  poly_scale_montgomery(&g,&g,4128752);
  poly_add(&g,&g,&h);
  poly_sub(&f,&f,&g);
  poly_freeze(&f);
  for(i=0;i<N;i++)
    printf("%d, ",f.coeffs[i]);
  printf("\n\n");

  return 0;
}
