#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include "randombytes.h"
#include "aes256ctr.h"
#include "comm.h"
#include "params.h"
#include "polyvec.h"
#include "poly.h"
#include "consts.h"

void expand_commkey(commkey *ck, const uint8_t rho[SYMBYTES]) {
  int i,j;
  aes256ctr_ctx state;

  aes256ctr_init(&state,rho,0);
  for(i=0;i<K;i++)
    for(j=0;j<L;j++) {
      aes256ctr_select(&state,(i << 16) + j);
      poly_uniform_preinit(&ck->b0[i].vec[j], &state);
    }

  for(i=0;i<K;i++)
    for(j=0;j<M;j++) {
      aes256ctr_select(&state,((K+i) << 16) + j);
      poly_uniform_preinit(&ck->bt[i].vec[j], &state);
    }

  for(i=0;i<M;i++)
    for(j=0;j<L;j++) {
      aes256ctr_select(&state,((2*K+i) << 16) + j);
      poly_uniform_preinit(&ck->bm[i].vec[j], &state);
    }
}

void commit(comm *t, commrnd *r, const polyvecm *msg, const commkey *ck) {
  int i;
  uint64_t nonce = 0;
  uint8_t seed[SYMBYTES];
  polyveck tag;
  aes256ctr_ctx state;

  randombytes(seed,SYMBYTES);
  aes256ctr_init(&state,seed,0);
  for(i=0;i<L;i++) {
    aes256ctr_select(&state,nonce++);
    poly_trinary_preinit(&r->s.vec[i],&state);
  }
  for(i=0;i<K;i++) {
    aes256ctr_select(&state,nonce++);
    poly_trinary_preinit(&r->e.vec[i],&state);
  }
  for(i=0;i<M;i++) {
    aes256ctr_select(&state,nonce++);
    poly_trinary_preinit(&r->em.vec[i],&state);
  }

  polyvecl_ntt(&r->s);
  polyveck_ntt(&r->e);
  polyvecm_ntt(&r->em);

  for(i=0;i<K;i++)
    polyvecl_pointwise_acc_montgomery(&t->t0.vec[i],&ck->b0[i],&r->s);
  for(i=0;i<K;i++)
    polyvecm_pointwise_acc_montgomery(&tag.vec[i],&ck->bt[i],&r->em);
  for(i=0;i<M;i++)
    polyvecl_pointwise_acc_montgomery(&t->tm.vec[i],&ck->bm[i],&r->s);

  polyveck_add(&t->t0,&t->t0,&tag);
  polyveck_scale_montgomery(&t->t0,&t->t0,MONTSQ);
  polyveck_add(&t->t0,&t->t0,&r->e);
  polyvecm_scale_montgomery(&t->tm,&t->tm,MONTSQ);
  polyvecm_add(&t->tm,&t->tm,&r->em);
  polyvecm_add(&t->tm,&t->tm,msg);
}
