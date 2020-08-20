#include <stdint.h>
#include <string.h>
#include "aes256ctr.h"
#include "params.h"
#include "consts.h"
#include "comm.h"
#include "apprshort.h"
#include "polyvec.h"
#include "poly.h"
#include "linear.h"

extern const polyvecm g[R];
static poly apprshort_gamma[R];
static poly apprshort_atgamma[R];

static void poly_invntt_transpose(poly *r, const poly *a) {
  int i;

  r->coeffs[0] = a->coeffs[0];
  for(i=1;i<N;i++)
    r->coeffs[i] = -a->coeffs[N-i];
  poly_ntt(r);
  poly_scale_montgomery(r,r,33554432);
}

static inline void poly_shift(poly *r, const poly *a, int i) {
  if(i == 0) *r = *a;
  else if(i == 1) poly_pointwise_montgomery(r,a,&nttx);
  else if(i == 2) poly_pointwise_montgomery(r,a,&nttx2);
  else if(i == 3) poly_pointwise_montgomery(r,a,&nttx3);
}

void linear_challenge_prehash(const uint8_t chash[SYMBYTES]) {
  int i;
  uint64_t nonce = 0;
  aes256ctr_ctx state;

  aes256ctr_init(&state,chash,0);
  for(i=0;i<R;i++) {
    aes256ctr_select(&state,nonce++);
    poly_uniform_preinit(&apprshort_gamma[i],&state);
    apprshort_challenge_mul(&apprshort_atgamma[i],&apprshort_gamma[i]);
  }
}

void linear(poly *h, poly vprime[R], const polyvecm *msg, const poly *l, uint8_t chash[SYMBYTES]) {
  int i,j,k;
  uint64_t nonce = 0;
  poly gamma, atgamma, btgamma;
  poly tmp, tmp2;
  aes256ctr_ctx aesctx;

  aes256ctr_init(&aesctx,chash,0);
  *h = msg->vec[M-1];
  for(i=0;i<R;i++)
    vprime[i] = g[i].vec[M-1];

  /* a = V^-1\hat{a}, b = V^-1\hat{b} */
  for(i=0;i<2;i++) {
    for(j=0;j<R;j++) {
      aes256ctr_select(&aesctx,nonce++);
      poly_uniform_preinit(&gamma,&aesctx);
      poly_invntt_transpose(&atgamma,&gamma);

      for(k=0;k<R;k++) {
        poly_pointwise_montgomery(&tmp,&gamma,&g[k].vec[i]);
        poly_pointwise_montgomery(&tmp2,&atgamma,&g[k].vec[3+i]);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_trace65_ntt(&tmp,&tmp);
        poly_shift(&tmp,&tmp,j);
        poly_add(&vprime[k],&vprime[k],&tmp);
      }

      poly_pointwise_montgomery(&tmp,&gamma,&msg->vec[i]);
      poly_pointwise_montgomery(&tmp2,&atgamma,&msg->vec[3+i]);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,j);
      poly_add(h,h,&tmp);
    }
  }

  /* c - V^-1\hat{ab} = (X-2)f */
  for(i=0;i<R;i++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&gamma,&aesctx);
    poly_invntt_transpose(&atgamma,&gamma);
    for(j=0;j<N-1;j++)
      btgamma.coeffs[j] = -2*gamma.coeffs[j] + gamma.coeffs[j+1];
    btgamma.coeffs[N-1] = -2*gamma.coeffs[N-1];

    for(j=0;j<R;j++) {
      poly_pointwise_montgomery(&tmp,&gamma,&g[j].vec[2]);
      poly_pointwise_montgomery(&tmp2,&atgamma,&g[j].vec[5]);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_pointwise_montgomery(&tmp2,&btgamma,&g[j].vec[6]);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,i);
      poly_add(&vprime[j],&vprime[j],&tmp);
    }

    poly_pointwise_montgomery(&tmp,&gamma,&msg->vec[2]);
    poly_pointwise_montgomery(&tmp2,&atgamma,&msg->vec[5]);
    poly_sub(&tmp,&tmp,&tmp2);
    poly_pointwise_montgomery(&tmp2,&btgamma,&msg->vec[6]);
    poly_sub(&tmp,&tmp,&tmp2);
    poly_trace65_ntt(&tmp,&tmp);
    poly_shift(&tmp,&tmp,i);
    poly_add(h,h,&tmp);
  }

  /* l = y_l + Bf */
  for(i=0;i<R;i++) {
    for(j=0;j<R;j++) {
      poly_pointwise_montgomery(&tmp,&apprshort_atgamma[i],&g[j].vec[6]);
      poly_pointwise_montgomery(&tmp2,&apprshort_gamma[i],&g[j].vec[M-3]);
      poly_add(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,i);
      poly_add(&vprime[j],&vprime[j],&tmp);
    }

    poly_pointwise_montgomery(&tmp,&apprshort_atgamma[i],&msg->vec[6]);
    poly_sub(&tmp2,&msg->vec[M-3],l);
    poly_pointwise_montgomery(&tmp2,&apprshort_gamma[i],&tmp2);
    poly_add(&tmp,&tmp,&tmp2);
    poly_trace65_ntt(&tmp,&tmp);
    poly_shift(&tmp,&tmp,i);
    poly_add(h,h,&tmp);
  }

  poly_invntt(h);
  poly_reduce(h);
  poly_freeze(h);
  for(i=0;i<R;i++)
    poly_freeze(&vprime[i]);
}

int linear_verify(poly vprime[R], const uint8_t chash[SYMBYTES], const poly *h, const poly *l,
                  const poly c[R], const commrnd z[R], const comm *tp, const commkey *ckp)
{
  int i,j,k;
  uint64_t nonce = 0;
  polyvecl zshat[R];
  poly chat[R], hhat;
  poly gamma[R], atgamma[R], btgamma[R];
  poly vtmp1, vtmp2, vtmp3, tmp, tmp2;
  aes256ctr_ctx aesctx;

  if(h->coeffs[0] || h->coeffs[1] || h->coeffs[2] || h->coeffs[3])
    return 1;

  for(j=0;j<R;j++) {
    zshat[j] = z[j].s;
    chat[j] = c[j];
    polyvecl_ntt(&zshat[j]);
    poly_ntt(&chat[j]);
  }

  aes256ctr_init(&aesctx,chash,0);
  hhat = *h;
  poly_ntt(&hhat);
  for(j=0;j<R;j++) {
    polyvecl_pointwise_acc_montgomery(&vprime[j],&ckp->bm[M-1],&zshat[j]);
    poly_sub(&tmp,&tp->tm.vec[M-1],&hhat);
    poly_pointwise_montgomery(&tmp,&chat[j],&tmp);
    poly_sub(&vprime[j],&vprime[j],&tmp);
    poly_scale_montgomery(&vprime[j],&vprime[j],MONTSQ);
    tmp = z[j].em.vec[M-1];
    poly_ntt(&tmp);
    poly_add(&vprime[j],&vprime[j],&tmp);
  }

  for(i=0;i<2;i++) {
    for(k=0;k<R;k++) {
      aes256ctr_select(&aesctx,nonce++);
      poly_uniform_preinit(&gamma[k],&aesctx);
      poly_invntt_transpose(&atgamma[k],&gamma[k]);
    }

    for(j=0;j<R;j++) {
      polyvecl_pointwise_acc_montgomery(&vtmp1,&ckp->bm[i],&zshat[j]);
      polyvecl_pointwise_acc_montgomery(&vtmp2,&ckp->bm[3+i],&zshat[j]);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[i]);
      poly_sub(&vtmp1,&vtmp1,&tmp);
      poly_scale_montgomery(&vtmp1,&vtmp1,MONTSQ);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[3+i]);
      poly_sub(&vtmp2,&vtmp2,&tmp);
      poly_scale_montgomery(&vtmp2,&vtmp2,MONTSQ);
      tmp = z[j].em.vec[i];
      poly_ntt(&tmp);
      poly_add(&vtmp1,&vtmp1,&tmp);
      tmp = z[j].em.vec[3+i];
      poly_ntt(&tmp);
      poly_add(&vtmp2,&vtmp2,&tmp);

      for(k=0;k<R;k++) {
        poly_pointwise_montgomery(&tmp,&gamma[k],&vtmp1);
        poly_pointwise_montgomery(&tmp2,&atgamma[k],&vtmp2);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_trace65_ntt(&tmp,&tmp);
        poly_shift(&tmp,&tmp,k);
        poly_add(&vprime[j],&vprime[j],&tmp);
      }
    }
  }

  for(k=0;k<R;k++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&gamma[k],&aesctx);
    poly_invntt_transpose(&atgamma[k],&gamma[k]);
    for(j=0;j<N-1;j++)
      btgamma[k].coeffs[j] = -2*gamma[k].coeffs[j] + gamma[k].coeffs[j+1];
    btgamma[k].coeffs[N-1] = -2*gamma[k].coeffs[j];
  }

  for(j=0;j<R;j++) {
    polyvecl_pointwise_acc_montgomery(&vtmp1,&ckp->bm[2],&zshat[j]);
    polyvecl_pointwise_acc_montgomery(&vtmp2,&ckp->bm[5],&zshat[j]);
    polyvecl_pointwise_acc_montgomery(&vtmp3,&ckp->bm[6],&zshat[j]);
    poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[2]);
    poly_sub(&vtmp1,&vtmp1,&tmp);
    poly_scale_montgomery(&vtmp1,&vtmp1,MONTSQ);
    poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[5]);
    poly_sub(&vtmp2,&vtmp2,&tmp);
    poly_scale_montgomery(&vtmp2,&vtmp2,MONTSQ);
    poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[6]);
    poly_sub(&vtmp3,&vtmp3,&tmp);
    poly_scale_montgomery(&vtmp3,&vtmp3,MONTSQ);
    tmp = z[j].em.vec[2];
    poly_ntt(&tmp);
    poly_add(&vtmp1,&vtmp1,&tmp);
    tmp = z[j].em.vec[5];
    poly_ntt(&tmp);
    poly_add(&vtmp2,&vtmp2,&tmp);
    tmp = z[j].em.vec[6];
    poly_ntt(&tmp);
    poly_add(&vtmp3,&vtmp3,&tmp);

    for(k=0;k<R;k++) {
      poly_pointwise_montgomery(&tmp,&gamma[k],&vtmp1);
      poly_pointwise_montgomery(&tmp2,&atgamma[k],&vtmp2);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_pointwise_montgomery(&tmp2,&btgamma[k],&vtmp3);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,k);
      poly_add(&vprime[j],&vprime[j],&tmp);
    }
  }

  for(j=0;j<R;j++) {
    polyvecl_pointwise_acc_montgomery(&vtmp1,&ckp->bm[6],&zshat[j]);
    polyvecl_pointwise_acc_montgomery(&vtmp2,&ckp->bm[M-3],&zshat[j]);
    poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[6]);
    poly_sub(&vtmp1,&vtmp1,&tmp);
    poly_scale_montgomery(&vtmp1,&vtmp1,MONTSQ);
    poly_sub(&tmp,&tp->tm.vec[M-3],l);
    poly_pointwise_montgomery(&tmp,&chat[j],&tmp);
    poly_sub(&vtmp2,&vtmp2,&tmp);
    poly_scale_montgomery(&vtmp2,&vtmp2,MONTSQ);
    tmp = z[j].em.vec[6];
    poly_ntt(&tmp);
    poly_add(&vtmp1,&vtmp1,&tmp);
    tmp = z[j].em.vec[M-3];
    poly_ntt(&tmp);
    poly_add(&vtmp2,&vtmp2,&tmp);

    for(k=0;k<R;k++) {
      poly_pointwise_montgomery(&tmp,&apprshort_atgamma[k],&vtmp1);
      poly_pointwise_montgomery(&tmp2,&apprshort_gamma[k],&vtmp2);
      poly_add(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,k);
      poly_add(&vprime[j],&vprime[j],&tmp);
    }
  }

  for(j=0;j<R;j++)
    poly_freeze(&vprime[j]);

  return 0;
}
