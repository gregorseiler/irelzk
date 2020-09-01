#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include "fips202.h"
#include "randombytes.h"
#include "params.h"
#include "consts.h"
#include "main.h"
#include "opening.h"
#include "product.h"
#include "linear.h"
#include "comm.h"
#include "polyvec.h"
#include "poly.h"

void main_proof(proof *p, comm *t, commrnd *r, const uint8_t rho[SYMBYTES],
                poly At[8][8], const poly u[8], const poly s[16])
{
  int i, rej;
  uint8_t symbuf[2*SYMBYTES+SHAKE128_RATE];
  uint8_t *seed = symbuf;
  uint8_t *thash = symbuf+SYMBYTES;
  uint8_t *chash = symbuf+2*SYMBYTES;
  uint16_t nonce = 0;
  commkey ck;
  polyveck t0low, w1[R];
  polyvecm msg;
  poly v, vprime[R], tmp, tmp2;
  keccak_state kecst;

  randombytes(seed,SYMBYTES);

  for(i=0;i<16;i++)
    msg.vec[i] = s[i];
  memset(&msg.vec[M-3],0,sizeof(poly));
  memset(&msg.vec[M-2],0,sizeof(poly));
  linear_first(&msg.vec[M-1]);

  expand_commkey(&ck,rho);
  commit(t,r,&msg,&ck);
  polyveck_invntt(&t->t0);
  polyveck_power2round(&t->t0,&t0low,&t->t0);
  polyveck_ntt(&t->t0);
  polyveck_ntt(&t0low);

  shake128_init(&kecst);
  shake128_absorb(&kecst,rho,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->t0,sizeof(polyveck));
  shake128_absorb(&kecst,(uint8_t *)&t->tm,(M-3)*sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-1],sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  linear_challenge_prehash(thash,At);
  linear_last(&p->h,&msg,u);

  shake128_init(&kecst);
  shake128_absorb(&kecst,thash,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&p->h,sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  opening_init(&t0low,r,&ck);
  do {
    opening_first(w1,seed,nonce);
    nonce += R*(L+K+M);

    shake128_init(&kecst);
    shake128_absorb(&kecst,thash,SYMBYTES);
    shake128_absorb(&kecst,(uint8_t *)w1,R*sizeof(polyveck));
    shake128_finalize(&kecst);
    shake128_squeezeblocks(chash,1,&kecst);

    product(&v,&msg,chash);
    poly_add(&tmp,&t->tm.vec[M-3],&msg.vec[M-3]);
    poly_add(&tmp2,&t->tm.vec[M-2],&msg.vec[M-2]);
    linear(vprime);

    shake128_init(&kecst);
    shake128_absorb(&kecst,chash,SYMBYTES);
    shake128_absorb(&kecst,(uint8_t *)&tmp,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)&tmp2,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)&v,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)vprime,R*sizeof(poly));
    shake128_finalize(&kecst);
    shake128_squeezeblocks(chash,1,&kecst);

    challenge_prehash(p->c,chash);
    rej = opening_last(p->z,p->c,w1);
  } while(rej);

  t->tm.vec[M-3] = tmp;
  t->tm.vec[M-2] = tmp2;
}

int main_proof_verify(const proof *p, const comm *t, const uint8_t rho[SYMBYTES],
                      poly At[8][8], const poly u[8])
{
  uint8_t symbuf[SYMBYTES+SHAKE128_RATE];
  uint8_t *thash = symbuf;
  uint8_t *chash = symbuf+SYMBYTES;
  commkey ck;
  poly v, vprime[R];
  polyveck w1[R];
  keccak_state kecst;

  shake128_init(&kecst);
  shake128_absorb(&kecst,rho,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->t0,sizeof(polyveck));
  shake128_absorb(&kecst,(uint8_t *)&t->tm,(M-3)*sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-1],sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  linear_challenge_prehash(thash,At);
  if(linear_verify_first(&p->h))
    return 1;

  shake128_init(&kecst);
  shake128_absorb(&kecst,thash,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&p->h,sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  expand_commkey(&ck,rho);
  if(opening_verify_first(w1,p->c,p->z,t,&ck))
    return 1;

  shake128_init(&kecst);
  shake128_absorb(&kecst,thash,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)w1,R*sizeof(polyveck));
  shake128_finalize(&kecst);
  shake128_squeezeblocks(chash,1,&kecst);

  if(product_verify(&v,chash,p->c,p->z,t,&ck))
    return 1;
  if(linear_verify(vprime,&p->h,u,p->c,p->z,t,&ck))
    return 1;

  shake128_init(&kecst);
  shake128_absorb(&kecst,chash,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-3],sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-2],sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&v,sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)vprime,4*sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeezeblocks(chash,1,&kecst);

  if(opening_verify_last(p->c,chash))
    return 1;

  return 0;
}
