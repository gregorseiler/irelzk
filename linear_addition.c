#include <stdint.h>
#include <string.h>
#include "aes256ctr.h"
#include "params.h"
#include "consts.h"
#include "comm.h"
#include "polyvec.h"
#include "poly.h"
#include "linear.h"

extern const polyvecm g[R];

static inline void poly_shift(poly *r, const poly *a, int i) {
  if(i == 0) *r = *a;
  else if(i == 1) poly_pointwise_montgomery(r,a,&nttx);
  else if(i == 2) poly_pointwise_montgomery(r,a,&nttx2);
  else if(i == 3) poly_pointwise_montgomery(r,a,&nttx3);
}

void linear(poly *h, poly vprime[R], const polyvecm *msg, uint8_t chash[SYMBYTES]) {
  int i,j;
  uint64_t nonce = 0;
  poly gamma, atgamma;
  poly vtmp[R], mpr, tmp, tmp2;
  aes256ctr_ctx aesctx;

  for(j=0;j<R;j++) {
    poly_add(&vtmp[j],&g[j].vec[0],&g[j].vec[1]);
    poly_sub(&vtmp[j],&vtmp[j],&g[j].vec[2]);
  }

  poly_add(&mpr,&msg->vec[0],&msg->vec[1]);
  poly_sub(&mpr,&mpr,&msg->vec[2]);

  aes256ctr_init(&aesctx,chash,0);
  *h = msg->vec[M-1];
  for(j=0;j<R;j++)
    vprime[j] = g[j].vec[M-1];

  for(i=0;i<R;i++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&gamma,&aesctx);
    for(j=0;j<N-1;j++)
      atgamma.coeffs[j] = 2*gamma.coeffs[j] - gamma.coeffs[j+1];
    atgamma.coeffs[N-1] = 2*gamma.coeffs[N-1];

    for(j=0;j<4;j++) {
      poly_pointwise_montgomery(&tmp,&gamma,&vtmp[j]);
      poly_pointwise_montgomery(&tmp2,&atgamma,&g[j].vec[3]);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,i);
      poly_add(&vprime[j],&vprime[j],&tmp);
    }

    poly_pointwise_montgomery(&tmp,&gamma,&mpr);
    poly_pointwise_montgomery(&tmp2,&atgamma,&msg->vec[3]);
    poly_sub(&tmp,&tmp,&tmp2);
    poly_trace65_ntt(&tmp,&tmp);
    poly_shift(&tmp,&tmp,i);
    poly_add(h,h,&tmp);
  }

  poly_invntt(h);
  poly_reduce(h);
  poly_freeze(h);
  for(j=0;j<R;j++)
    poly_freeze(&vprime[j]);
}

int linear_verify(poly vprime[R], const uint8_t chash[SYMBYTES], const poly *h,
                  const poly c[R], const commrnd z[R], const comm *tp, const commkey *ckp)
{
  int i,j;
  uint64_t nonce = 0;
  poly gamma[R], atgamma[R];
  polyvecl bpr, zshat;
  poly tpr, vtmp1, vtmp2, tmp, chat, hhat;
  aes256ctr_ctx aesctx;

  if(h->coeffs[0] || h->coeffs[1] || h->coeffs[2] || h->coeffs[3])
    return 1;

  aes256ctr_init(&aesctx,chash,0);
  for(i=0;i<R;i++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&gamma[i],&aesctx);
    for(j=0;j<N-1;j++)
      atgamma[i].coeffs[j] = 2*gamma[i].coeffs[j] - gamma[i].coeffs[j+1];
    atgamma[i].coeffs[N-1] = 2*gamma[i].coeffs[N-1];
  }

  polyvecl_add(&bpr,&ckp->bm[0],&ckp->bm[1]);
  polyvecl_sub(&bpr,&bpr,&ckp->bm[2]);

  poly_add(&tpr,&tp->tm.vec[0],&tp->tm.vec[1]);
  poly_sub(&tpr,&tpr,&tp->tm.vec[2]);

  hhat = *h;
  poly_ntt(&hhat);

  for(i=0;i<4;i++) {
    zshat = z[i].s;
    polyvecl_ntt(&zshat);
    polyvecl_pointwise_acc_montgomery(&vtmp1,&bpr,&zshat);
    polyvecl_pointwise_acc_montgomery(&vtmp2,&ckp->bm[3],&zshat);
    polyvecl_pointwise_acc_montgomery(&vprime[i],&ckp->bm[M-1],&zshat);

    chat = c[i];
    poly_ntt(&chat);
    poly_pointwise_montgomery(&tmp,&chat,&tpr);
    poly_sub(&vtmp1,&vtmp1,&tmp);
    poly_scale_montgomery(&vtmp1,&vtmp1,MONTSQ);
    poly_pointwise_montgomery(&tmp,&chat,&tp->tm.vec[3]);
    poly_sub(&vtmp2,&vtmp2,&tmp);
    poly_scale_montgomery(&vtmp2,&vtmp2,MONTSQ);
    poly_sub(&tmp,&tp->tm.vec[5],&hhat);
    poly_pointwise_montgomery(&tmp,&chat,&tmp);
    poly_sub(&vprime[i],&vprime[i],&tmp);
    poly_scale_montgomery(&vprime[i],&vprime[i],MONTSQ);

    poly_add(&tmp,&z[i].em.vec[0],&z[i].em.vec[1]);
    poly_sub(&tmp,&tmp,&z[i].em.vec[2]);
    poly_ntt(&tmp);
    poly_add(&vtmp1,&vtmp1,&tmp);

    tmp = z[i].em.vec[3];
    poly_ntt(&tmp);
    poly_add(&vtmp2,&vtmp2,&tmp);

    tmp = z[i].em.vec[5];
    poly_ntt(&tmp);
    poly_add(&vprime[i],&vprime[i],&tmp);

    for(j=0;j<R;j++) {
      poly_pointwise_montgomery(&tmp,&gamma[j],&vtmp1);
      poly_pointwise_montgomery(&chat,&atgamma[j],&vtmp2);
      poly_sub(&tmp,&tmp,&chat);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,j);
      poly_add(&vprime[i],&vprime[i],&tmp);
    }

    poly_freeze(&vprime[i]);
  }

  return 0;
}
