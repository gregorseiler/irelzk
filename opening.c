#include <stdint.h>
#include <string.h>
#include "fips202.h"
#include "randombytes.h"
#include "params.h"
#include "opening.h"
#include "comm.h"
#include "polyvec.h"
#include "consts.h"

static const polyveck *t0low;
static const commrnd *r;
static const commkey *ck;
static commrnd y[R];
polyvecm g[R];

void challenge_prehash(poly c[R], const uint8_t chash[N/4]) {
  int i, j;
  uint8_t b;
  const int32_t lut[4] = {0, 0, 1, -1};

  memset(c,0,R*sizeof(poly));
  for(i=0;i<R;i++) {
    for(j=0;j<N/(4*R);j++) {
      b = chash[N/(4*R)*i+j];
      c[i].coeffs[16*j+0*R] = lut[(b >> 0)&3];
      c[i].coeffs[16*j+1*R] = lut[(b >> 2)&3];
      c[i].coeffs[16*j+2*R] = lut[(b >> 4)&3];
      c[i].coeffs[16*j+3*R] = lut[(b >> 6)&3];
    }
  }
}

void challenge(poly c[R], const polyveck w[R]) {
  uint8_t chash[N/4];
  shake128(chash,N/4,(uint8_t *)w,R*sizeof(polyveck));
  challenge_prehash(c,chash);
}

void opening_init(const polyveck *t0lowp, const commrnd *rp, const commkey *ckp) {
  t0low = t0lowp;
  r = rp;
  ck = ckp;
}

void opening_first(polyveck w1[4], const uint8_t seed[SYMBYTES], uint16_t nonce) {
  int i,j;
  poly tmp;
  aes256ctr_ctx state;

  aes256ctr_init(&state,seed,nonce);
  for(i=0;i<R;i++) {
    for(j=0;j<L;j++) {
      aes256ctr_select(&state,nonce++);
      poly_uniform_gamma_preinit(&y[i].s.vec[j],&state);
    }
    for(j=0;j<M;j++) {
      aes256ctr_select(&state,nonce++);
      poly_uniform_gamma_preinit(&y[i].em.vec[j],&state);
    }

    polyvecl_ntt(&y[i].s);
    polyvecm_ntt(&y[i].em);

    for(j=0;j<K;j++)
      polyvecl_pointwise_acc_montgomery(&w1[i].vec[j],&ck->b0[j],&y[i].s);
    for(j=0;j<K;j++) {
      polyvecm_pointwise_acc_montgomery(&tmp,&ck->bt[j],&y[i].em);
      poly_add(&w1[i].vec[j],&w1[i].vec[j],&tmp);
    }
    polyveck_invntt_tomont(&w1[i]);
    polyveck_decompose(&w1[i],&y[i].e,&w1[i]);

    for(j=0;j<M;j++)
      polyvecl_pointwise_acc_montgomery(&g[i].vec[j],&ck->bm[j],&y[i].s);
    polyvecm_scale_montgomery(&g[i],&g[i],MONTSQ);
    polyvecm_add(&g[i],&g[i],&y[i].em);
  }
}

int opening_last(commrnd z[R], const poly c[R], const polyveck w1[R]) {
  int i,j;
  poly chat[R];

  for(i=0;i<R;i++) {
    chat[i] = c[i];
    poly_ntt(&chat[i]);
    poly_scale_montgomery(&chat[i],&chat[i],MONTSQ);
  }

  for(i=0;i<R;i++) {
    for(j=0;j<L;j++)
      poly_pointwise_montgomery(&z[i].s.vec[j],&chat[i],&r->s.vec[j]);
    polyvecl_add(&z[i].s,&z[i].s,&y[i].s);
    polyvecl_invntt(&z[i].s);
    for(j=0;j<L;j++)
      poly_reduce(&z[i].s.vec[j]);
    if(polyvecl_chknorm(&z[i].s,GAMMA1-BETA))
      return 1;

    for(j=0;j<M;j++)
      poly_pointwise_montgomery(&z[i].em.vec[j],&chat[i],&r->em.vec[j]);
    polyvecm_add(&z[i].em,&z[i].em,&y[i].em);
    polyvecm_invntt(&z[i].em);
    for(j=0;j<M;j++)
      poly_reduce(&z[i].em.vec[j]);
    if(polyvecm_chknorm(&z[i].em,GAMMA1-BETA))
      return 1;

    for(j=0;j<K;j++)
      poly_pointwise_montgomery(&z[i].e.vec[j],&chat[i],&r->e.vec[j]);
    polyveck_invntt(&z[i].e);
    polyveck_sub(&y[i].e,&y[i].e,&z[i].e);
    if(polyveck_chknorm(&y[i].e,GAMMA2-BETA))
      return 1;
  }

  for(i=0;i<R;i++) {
    for(j=0;j<K;j++)
      poly_pointwise_montgomery(&z[i].e.vec[j],&chat[i],&t0low->vec[j]);
    polyveck_invntt(&z[i].e);
    if(polyveck_chknorm(&z[i].e,GAMMA2)) return 1;
    polyveck_add(&y[i].e,&y[i].e,&z[i].e);
    polyveck_makehint(&z[i].e,&w1[i],&y[i].e);
  }

  return 0;
}

int opening_verify_first(polyveck w1[R], const poly c[R], const commrnd z[R],
                         const comm *tp, const commkey *ckp)
{
  int i,j;
  polyvecl zshat;
  polyvecm zmhat;
  poly tmp, chat;

  for(i=0;i<R;i++) {
    if(polyvecl_chknorm(&z[i].s,GAMMA1-BETA))
      return 1;
    if(polyvecm_chknorm(&z[i].em,GAMMA1-BETA))
      return 1;
  }

  for(i=0;i<R;i++) {
    zshat = z[i].s;
    zmhat = z[i].em;
    polyvecl_ntt(&zshat);
    polyvecm_ntt(&zmhat);
    for(j=0;j<K;j++)
      polyvecl_pointwise_acc_montgomery(&w1[i].vec[j],&ckp->b0[j],&zshat);
    for(j=0;j<K;j++) {
      polyvecm_pointwise_acc_montgomery(&tmp,&ckp->bt[j],&zmhat);
      poly_add(&w1[i].vec[j],&w1[i].vec[j],&tmp);
    }

    chat = c[i];
    poly_ntt(&chat);
    for(j=0;j<K;j++) {
      poly_scale_montgomery(&tmp,&tp->t0.vec[j],4128752);
      poly_pointwise_montgomery(&tmp,&chat,&tmp);
      poly_sub(&w1[i].vec[j],&w1[i].vec[j],&tmp);
    }

    polyveck_invntt_tomont(&w1[i]);
    polyveck_usehint(&w1[i],&w1[i],&z[i].e);
  }

  return 0;
}

int opening_verify_last(const poly c[R], const uint8_t chash[N/4]) {
  int i,j;
  poly c2[R];

  challenge_prehash(c2,chash);
  for(i=0;i<R;i++)
    for(j=0;j<N;j++)
      if(c[i].coeffs[j] != c2[i].coeffs[j])
        return 1;

  return 0;
}
