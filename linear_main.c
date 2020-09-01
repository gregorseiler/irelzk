#include <stdint.h>
#include <string.h>
#include "aes256ctr.h"
#include "randombytes.h"
#include "params.h"
#include "consts.h"
#include "comm.h"
#include "apprshort.h"
#include "polyvec.h"
#include "poly.h"
#include "linear.h"

extern const polyvecm g[R];
static poly gammavec[R][8];
static poly atgammavec[R][8];

static inline void poly_shift(poly *r, const poly *a, int i) {
  if(i == 0) *r = *a;
  else if(i == 1) poly_pointwise_montgomery(r,a,&nttx);
  else if(i == 2) poly_pointwise_montgomery(r,a,&nttx2);
  else if(i == 3) poly_pointwise_montgomery(r,a,&nttx3);
}

void linear_first(poly *gp) {
  int i;
  uint8_t seed[SYMBYTES];

  randombytes(seed,SYMBYTES);
  poly_uniform(gp,seed,0);
  for(i=0;i<R;i++)
    gp->coeffs[i] = 0;
  poly_ntt(gp);
}

void linear_challenge_prehash(const uint8_t chash[SYMBYTES], poly At[8][8]) {
  int i,j,k;
  uint64_t nonce = 0;
  poly ghat[8], tmp;
  aes256ctr_ctx state;

  aes256ctr_init(&state,chash,0);
  for(i=0;i<R;i++) {
    for(j=0;j<8;j++) {
      aes256ctr_select(&state,nonce++);
      poly_uniform_preinit(&gammavec[i][j],&state);
    }

    for(j=0;j<8;j++) {
      ghat[j] = gammavec[i][j];
      poly_ntt(&ghat[j]);
    }
    for(j=0;j<8;j++) {
      for(k=0;k<8;k++) {
        poly_pointwise_montgomery(&tmp,&At[j][k],&ghat[k]);
        if(k == 0) atgammavec[i][j] = tmp;
        else poly_add(&atgammavec[i][j],&atgammavec[i][j],&tmp);
      }
      poly_invntt_tomont(&atgammavec[i][j]);
    }
  }
}

void linear_last(poly *h, const polyvecm *msg, const poly u[8]) {
  int i,j;
  poly tmp, tmp2;

  *h = msg->vec[M-1];
  for(i=0;i<R;i++) {
    for(j=0;j<8;j++) {
      poly_pointwise_montgomery(&tmp2,&atgammavec[i][j],&msg->vec[j]);
      if(j == 0) tmp = tmp2;
      else poly_add(&tmp,&tmp,&tmp2);

      poly_sub(&tmp2,&msg->vec[8+j],&u[j]);
      poly_pointwise_montgomery(&tmp2,&gammavec[i][j],&tmp2);
      poly_add(&tmp,&tmp,&tmp2);
    }
    poly_trace65_ntt(&tmp,&tmp);
    poly_shift(&tmp,&tmp,i);
    poly_add(h,h,&tmp);
  }
  poly_invntt(h);
  poly_reduce(h);
  poly_freeze(h);
}

void linear(poly vprime[R]) {
  int i,j,k;
  poly tmp, tmp2;

  for(i=0;i<R;i++) {
    vprime[i] = g[i].vec[M-1];
    for(j=0;j<R;j++) {
      for(k=0;k<8;k++) {
        poly_pointwise_montgomery(&tmp,&atgammavec[j][k],&g[i].vec[k]);
        poly_pointwise_montgomery(&tmp2,&gammavec[j][k],&g[i].vec[8+k]);
        poly_add(&tmp,&tmp,&tmp2);
        poly_trace65_ntt(&tmp,&tmp);
        poly_shift(&tmp,&tmp,j);
        poly_add(&vprime[i],&vprime[i],&tmp);
      }
    }
    poly_freeze(&vprime[i]);
  }
}

int linear_verify_first(const poly *h) {
  if(h->coeffs[0] || h->coeffs[1] || h->coeffs[2] || h->coeffs[3])
    return 1;

  return 0;
}

int linear_verify(poly vprime[R], const poly *h, const poly u[8], const poly c[R], const commrnd z[R],
                  const comm *tp, const commkey *ckp)
{
  int i,j,k;
  polyvecl zshat[R];
  poly chat[R], hhat;
  poly vtmp1, vtmp2, tmp, tmp2;

  for(i=0;i<R;i++) {
    zshat[i] = z[i].s;
    chat[i] = c[i];
    polyvecl_ntt(&zshat[i]);
    poly_ntt(&chat[i]);
  }

  hhat = *h;
  poly_ntt(&hhat);
  for(i=0;i<R;i++) {
    polyvecl_pointwise_acc_montgomery(&vprime[i],&ckp->bm[M-1],&zshat[i]);
    poly_sub(&tmp,&tp->tm.vec[M-1],&hhat);
    poly_pointwise_montgomery(&tmp,&chat[i],&tmp);
    poly_sub(&vprime[i],&vprime[i],&tmp);
    poly_scale_montgomery(&vprime[i],&vprime[i],MONTSQ);
    tmp = z[i].em.vec[M-1];
    poly_ntt(&tmp);
    poly_add(&vprime[i],&vprime[i],&tmp);
  }

  for(i=0;i<8;i++) {
    for(j=0;j<R;j++) {
      polyvecl_pointwise_acc_montgomery(&vtmp1,&ckp->bm[i],&zshat[j]);
      polyvecl_pointwise_acc_montgomery(&vtmp2,&ckp->bm[8+i],&zshat[j]);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[i]);
      poly_sub(&vtmp1,&vtmp1,&tmp);
      poly_scale_montgomery(&vtmp1,&vtmp1,MONTSQ);
      poly_sub(&tmp,&tp->tm.vec[8+i],&u[i]);
      poly_pointwise_montgomery(&tmp,&chat[j],&tmp);
      poly_sub(&vtmp2,&vtmp2,&tmp);
      poly_scale_montgomery(&vtmp2,&vtmp2,MONTSQ);
      tmp = z[j].em.vec[i];
      poly_ntt(&tmp);
      poly_add(&vtmp1,&vtmp1,&tmp);
      tmp = z[j].em.vec[8+i];
      poly_ntt(&tmp);
      poly_add(&vtmp2,&vtmp2,&tmp);

      for(k=0;k<R;k++) {
        poly_pointwise_montgomery(&tmp,&atgammavec[k][i],&vtmp1);
        poly_pointwise_montgomery(&tmp2,&gammavec[k][i],&vtmp2);
        poly_add(&tmp,&tmp,&tmp2);
        poly_trace65_ntt(&tmp,&tmp);
        poly_shift(&tmp,&tmp,k);
        poly_add(&vprime[j],&vprime[j],&tmp);
      }
    }
  }

  for(i=0;i<R;i++)
    poly_freeze(&vprime[i]);

  return 0;
}
