#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include "fips202.h"
#include "randombytes.h"
#include "params.h"
#include "consts.h"
#include "mult.h"
#include "opening.h"
#include "product.h"
#include "linear.h"
#include "apprshort.h"
#include "comm.h"
#include "polyvec.h"
#include "poly.h"

void multiplication_proof(proof *p, comm *t, commrnd *r, const uint8_t rho[SYMBYTES],
                          const uint64_t a[2], const uint64_t b[2])
{
  int i, rej;
  uint16_t nonce = 0;
  uint64_t c[2];
  uint8_t symbuf[2*SYMBYTES+SHAKE128_RATE];
  uint8_t *seed = symbuf;
  uint8_t *thash = symbuf+SYMBYTES;
  uint8_t *chash = symbuf+2*SYMBYTES;
  commkey ck;
  polyvecm msg;
  polyveck t0low, w1[R];
  poly v, vprime[R], tmp;
  keccak_state kecst;

  randombytes(seed,SYMBYTES);

  c[0] = _mulx_u64(a[0],b[0],(unsigned long long int *)&c[1]);
  memset(&msg,0,2*sizeof(poly));
  for(i=0;i<64;i++) {
    msg.vec[0].coeffs[i] = (a[0] >> i) & 1;
    msg.vec[1].coeffs[i] = (b[0] >> i) & 1;
    msg.vec[2].coeffs[i] = (c[0] >> i) & 1;
    msg.vec[2].coeffs[64+i] = (c[1] >> i) & 1;
  }

  msg.vec[3] = msg.vec[0];
  msg.vec[4] = msg.vec[1];
  poly_ntt(&msg.vec[3]);
  poly_ntt(&msg.vec[4]);
  poly_pointwise_montgomery(&msg.vec[6],&msg.vec[3],&msg.vec[4]);
  poly_scale_montgomery(&msg.vec[5],&msg.vec[6],MONTSQ);
  poly_invntt_tomont(&msg.vec[6]);
  poly_sub(&msg.vec[6],&msg.vec[2],&msg.vec[6]);
  poly_freeze(&msg.vec[6]);

  msg.vec[6].coeffs[0] = 0;
  for(i=1;i<128;i++)
    msg.vec[6].coeffs[i] = (msg.vec[6].coeffs[i-1] - msg.vec[6].coeffs[i]) >> 1;

  memset(&msg.vec[M-3],0,sizeof(poly));
  memset(&msg.vec[M-2],0,sizeof(poly));
  poly_uniform(&msg.vec[M-1],seed,nonce++);
  for(i=0;i<R;i++)
    msg.vec[M-1].coeffs[i] = 0;
  poly_ntt(&msg.vec[M-1]);

  expand_commkey(&ck,rho);
  commit(t,r,&msg,&ck);
  polyveck_invntt(&t->t0);
  polyveck_power2round(&t->t0,&t0low,&t->t0);
  polyveck_ntt(&t->t0);
  polyveck_ntt(&t0low);

  shake128_init(&kecst);
  shake128_absorb(&kecst,rho,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->t0,sizeof(polyveck));
  shake128_absorb(&kecst,(uint8_t *)&t->tm,7*sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-1],sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  do {
    apprshort_first(&msg.vec[M-3],seed,nonce++);
    poly_add(&tmp,&t->tm.vec[M-3],&msg.vec[M-3]);

    shake128_init(&kecst);
    shake128_absorb(&kecst,thash,SYMBYTES);
    shake128_absorb(&kecst,(uint8_t *)&tmp,sizeof(poly));
    shake128_finalize(&kecst);
    shake128_squeezeblocks(chash,1,&kecst);

    apprshort_challenge_prehash(chash);
    rej = apprshort_last(&p->l,&msg);
  } while(rej);

  t->tm.vec[M-3] = tmp;
  shake128_init(&kecst);
  shake128_absorb(&kecst,chash,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&p->l,sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  linear_challenge_prehash(thash);
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
    poly_add(&tmp,&t->tm.vec[M-2],&msg.vec[M-2]);
    linear(&p->h,vprime,&msg,&p->l,chash+SYMBYTES);

    shake128_init(&kecst);
    shake128_absorb(&kecst,chash,2*SYMBYTES);
    shake128_absorb(&kecst,(uint8_t *)&tmp,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)&v,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)&p->h,sizeof(poly));
    shake128_absorb(&kecst,(uint8_t *)vprime,R*sizeof(poly));
    shake128_finalize(&kecst);
    shake128_squeezeblocks(chash,1,&kecst);

    challenge_prehash(p->c,chash);
    rej = opening_last(p->z,p->c,w1);
  } while(rej);

  t->tm.vec[M-2] = tmp;
}

int multiplication_proof_verify(const proof *p, const comm *t, const uint8_t rho[SYMBYTES]) {
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
  shake128_absorb(&kecst,(uint8_t *)&t->tm,7*sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-1],sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  shake128_init(&kecst);
  shake128_absorb(&kecst,thash,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-3],sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeezeblocks(chash,1,&kecst);

  apprshort_challenge_prehash(chash);
  if(apprshort_verify(&p->l))
    return 1;

  shake128_init(&kecst);
  shake128_absorb(&kecst,chash,SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&p->l,sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeeze(thash,SYMBYTES,&kecst);

  linear_challenge_prehash(thash);
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
  if(linear_verify(vprime,chash+SYMBYTES,&p->h,&p->l,p->c,p->z,t,&ck))
    return 1;

  shake128_init(&kecst);
  shake128_absorb(&kecst,chash,2*SYMBYTES);
  shake128_absorb(&kecst,(uint8_t *)&t->tm.vec[M-2],sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&v,sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)&p->h,sizeof(poly));
  shake128_absorb(&kecst,(uint8_t *)vprime,4*sizeof(poly));
  shake128_finalize(&kecst);
  shake128_squeezeblocks(chash,1,&kecst);

  if(opening_verify_last(p->c,chash))
    return 1;

  return 0;
}