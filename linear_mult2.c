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

static void poly_invntt2_transpose(poly r[2], const poly a[2]) {
  int i;

  r[0].coeffs[0] = a[0].coeffs[0];
  r[1].coeffs[0] = a[1].coeffs[0];
  for(i=1;i<N;i++) {
    r[0].coeffs[i] = a[1].coeffs[N-i];
    r[1].coeffs[i] = a[0].coeffs[N-i];
  }
  poly_ntt2(r);
  poly_scale_montgomery(&r[0],&r[0],16777216);
  poly_scale_montgomery(&r[1],&r[1],16777216);
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

void linear(poly vprime[R], poly *h, const polyvecm *msg, const poly *l, uint8_t chash[SYMBYTES]) {
  int i,j,k;
  uint64_t nonce = 0;
  poly gamma[2], atgamma[2], btgamma[2];
  poly tmp, tmp2;
  aes256ctr_ctx aesctx;

  aes256ctr_init(&aesctx,chash,0);
  *h = msg->vec[M-1];
  for(i=0;i<R;i++)
    vprime[i] = g[i].vec[M-1];

  /* a = V^-1\hat{a}, b = V^-1\hat{b} */
  for(i=0;i<2;i++) {
    for(j=0;j<R;j++) {
      //FIXME: Only once
      aes256ctr_select(&aesctx,nonce++);
      poly_uniform_preinit(&gamma[0],&aesctx);
      aes256ctr_select(&aesctx,nonce++);
      poly_uniform_preinit(&gamma[1],&aesctx);
      poly_invntt2_transpose(atgamma,gamma);

      for(k=0;k<R;k++) {
        poly_pointwise_montgomery(&tmp,&gamma[0],&g[k].vec[i]);
        poly_pointwise_montgomery(&tmp2,&atgamma[0],&g[k].vec[4+2*i+0]);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_pointwise_montgomery(&tmp2,&atgamma[1],&g[k].vec[4+2*i+1]);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_trace65_ntt(&tmp,&tmp);
        poly_shift(&tmp,&tmp,j);
        poly_add(&vprime[k],&vprime[k],&tmp);
      }

      poly_pointwise_montgomery(&tmp,&gamma[0],&msg->vec[i]);
      poly_pointwise_montgomery(&tmp2,&atgamma[0],&msg->vec[4+2*i+0]);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_pointwise_montgomery(&tmp2,&atgamma[1],&msg->vec[4+2*i+1]);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,j);
      poly_add(h,h,&tmp);
    }
  }

  /* c - V^-1\hat{ab} = (X-2)f */
  for(i=0;i<R;i++) {
    //FIXME: Only once
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&gamma[0],&aesctx);
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&gamma[1],&aesctx);
    poly_invntt2_transpose(atgamma,gamma);
    for(j=0;j<N-1;j++) {
      btgamma[0].coeffs[j] = -2*gamma[0].coeffs[j] + gamma[0].coeffs[j+1];
      btgamma[1].coeffs[j] = -2*gamma[1].coeffs[j] + gamma[1].coeffs[j+1];
    }
    btgamma[0].coeffs[N-1] = -2*gamma[0].coeffs[N-1] + gamma[1].coeffs[0];
    btgamma[1].coeffs[N-1] = -2*gamma[1].coeffs[N-1];

    for(j=0;j<R;j++) {
      for(k=0;k<2;k++) {
        poly_pointwise_montgomery(&tmp,&gamma[k],&g[j].vec[2+k]);
        poly_pointwise_montgomery(&tmp2,&atgamma[k],&g[j].vec[8+k]);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_pointwise_montgomery(&tmp2,&btgamma[k],&g[j].vec[10+k]);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_trace65_ntt(&tmp,&tmp);
        poly_shift(&tmp,&tmp,i);
        poly_add(&vprime[j],&vprime[j],&tmp);
      }
    }

    for(k=0;k<2;k++) {
      poly_pointwise_montgomery(&tmp,&gamma[k],&msg->vec[2+k]);
      poly_pointwise_montgomery(&tmp2,&atgamma[k],&msg->vec[8+k]);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_pointwise_montgomery(&tmp2,&btgamma[k],&msg->vec[10+k]);
      poly_sub(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,i);
      poly_add(h,h,&tmp);
    }
  }

  /* l = y_l + Bf */
  for(i=0;i<R;i++) {
    for(j=0;j<R;j++) {
      poly_pointwise_montgomery(&tmp,&apprshort_atgamma[i],&g[j].vec[M-4]);
      poly_pointwise_montgomery(&tmp2,&apprshort_gamma[i],&g[j].vec[M-3]);
      poly_add(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,i);
      poly_add(&vprime[j],&vprime[j],&tmp);
    }

    poly_pointwise_montgomery(&tmp,&apprshort_atgamma[i],&msg->vec[M-4]);
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
  poly gamma[R][2], atgamma[R][2], btgamma[R][2];
  poly vtmp1, vtmp2, vtmp3, tmp, tmp2;
  aes256ctr_ctx aesctx;

  if(h->coeffs[0] || h->coeffs[1] || h->coeffs[2] || h->coeffs[3])
    return 1;

  for(i=0;i<R;i++) {
    zshat[i] = z[i].s;
    chat[i] = c[i];
    polyvecl_ntt(&zshat[i]);
    poly_ntt(&chat[i]);
  }

  aes256ctr_init(&aesctx,chash,0);
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

  for(i=0;i<2;i++) {
    for(j=0;j<R;j++) {
      aes256ctr_select(&aesctx,nonce++);
      poly_uniform_preinit(&gamma[j][0],&aesctx);
      aes256ctr_select(&aesctx,nonce++);
      poly_uniform_preinit(&gamma[j][1],&aesctx);
      poly_invntt2_transpose(atgamma[j],gamma[j]);
    }

    for(j=0;j<R;j++) {
      polyvecl_pointwise_acc_montgomery(&vtmp1,&ckp->bm[i],&zshat[j]);
      polyvecl_pointwise_acc_montgomery(&vtmp2,&ckp->bm[4+2*i+0],&zshat[j]);
      polyvecl_pointwise_acc_montgomery(&vtmp3,&ckp->bm[4+2*i+1],&zshat[j]);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[i]);
      poly_sub(&vtmp1,&vtmp1,&tmp);
      poly_scale_montgomery(&vtmp1,&vtmp1,MONTSQ);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[4+2*i+0]);
      poly_sub(&vtmp2,&vtmp2,&tmp);
      poly_scale_montgomery(&vtmp2,&vtmp2,MONTSQ);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[4+2*i+1]);
      poly_sub(&vtmp3,&vtmp3,&tmp);
      poly_scale_montgomery(&vtmp3,&vtmp3,MONTSQ);
      tmp = z[j].em.vec[i];
      poly_ntt(&tmp);
      poly_add(&vtmp1,&vtmp1,&tmp);
      tmp = z[j].em.vec[4+2*i+0];
      poly_ntt(&tmp);
      poly_add(&vtmp2,&vtmp2,&tmp);
      tmp = z[j].em.vec[4+2*i+1];
      poly_ntt(&tmp);
      poly_add(&vtmp3,&vtmp3,&tmp);

      for(k=0;k<R;k++) {
        poly_pointwise_montgomery(&tmp,&gamma[k][0],&vtmp1);
        poly_pointwise_montgomery(&tmp2,&atgamma[k][0],&vtmp2);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_pointwise_montgomery(&tmp2,&atgamma[k][1],&vtmp3);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_trace65_ntt(&tmp,&tmp);
        poly_shift(&tmp,&tmp,k);
        poly_add(&vprime[j],&vprime[j],&tmp);
      }
    }
  }

  for(i=0;i<R;i++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&gamma[i][0],&aesctx);
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&gamma[i][1],&aesctx);
    poly_invntt2_transpose(atgamma[i],gamma[i]);
    for(j=0;j<N-1;j++) {
      btgamma[i][0].coeffs[j] = -2*gamma[i][0].coeffs[j] + gamma[i][0].coeffs[j+1];
      btgamma[i][1].coeffs[j] = -2*gamma[i][1].coeffs[j] + gamma[i][1].coeffs[j+1];
    }
    btgamma[i][0].coeffs[N-1] = -2*gamma[i][0].coeffs[N-1] + gamma[i][1].coeffs[0];
    btgamma[i][1].coeffs[N-1] = -2*gamma[i][1].coeffs[N-1];
  }

  for(i=0;i<2;i++) {
    for(j=0;j<R;j++) {
      polyvecl_pointwise_acc_montgomery(&vtmp1,&ckp->bm[2+i],&zshat[j]);
      polyvecl_pointwise_acc_montgomery(&vtmp2,&ckp->bm[8+i],&zshat[j]);
      polyvecl_pointwise_acc_montgomery(&vtmp3,&ckp->bm[10+i],&zshat[j]);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[2+i]);
      poly_sub(&vtmp1,&vtmp1,&tmp);
      poly_scale_montgomery(&vtmp1,&vtmp1,MONTSQ);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[8+i]);
      poly_sub(&vtmp2,&vtmp2,&tmp);
      poly_scale_montgomery(&vtmp2,&vtmp2,MONTSQ);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[10+i]);
      poly_sub(&vtmp3,&vtmp3,&tmp);
      poly_scale_montgomery(&vtmp3,&vtmp3,MONTSQ);
      tmp = z[j].em.vec[2+i];
      poly_ntt(&tmp);
      poly_add(&vtmp1,&vtmp1,&tmp);
      tmp = z[j].em.vec[8+i];
      poly_ntt(&tmp);
      poly_add(&vtmp2,&vtmp2,&tmp);
      tmp = z[j].em.vec[10+i];
      poly_ntt(&tmp);
      poly_add(&vtmp3,&vtmp3,&tmp);

      for(k=0;k<R;k++) {
        poly_pointwise_montgomery(&tmp,&gamma[k][i],&vtmp1);
        poly_pointwise_montgomery(&tmp2,&atgamma[k][i],&vtmp2);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_pointwise_montgomery(&tmp2,&btgamma[k][i],&vtmp3);
        poly_sub(&tmp,&tmp,&tmp2);
        poly_trace65_ntt(&tmp,&tmp);
        poly_shift(&tmp,&tmp,k);
        poly_add(&vprime[j],&vprime[j],&tmp);
      }
    }
  }

  for(i=0;i<R;i++) {
    polyvecl_pointwise_acc_montgomery(&vtmp1,&ckp->bm[M-4],&zshat[i]);
    polyvecl_pointwise_acc_montgomery(&vtmp2,&ckp->bm[M-3],&zshat[i]);
    poly_pointwise_montgomery(&tmp,&chat[i],&tp->tm.vec[M-4]);
    poly_sub(&vtmp1,&vtmp1,&tmp);
    poly_scale_montgomery(&vtmp1,&vtmp1,MONTSQ);
    poly_sub(&tmp,&tp->tm.vec[M-3],l);
    poly_pointwise_montgomery(&tmp,&chat[i],&tmp);
    poly_sub(&vtmp2,&vtmp2,&tmp);
    poly_scale_montgomery(&vtmp2,&vtmp2,MONTSQ);
    tmp = z[i].em.vec[M-4];
    poly_ntt(&tmp);
    poly_add(&vtmp1,&vtmp1,&tmp);
    tmp = z[i].em.vec[M-3];
    poly_ntt(&tmp);
    poly_add(&vtmp2,&vtmp2,&tmp);

    for(j=0;j<R;j++) {
      poly_pointwise_montgomery(&tmp,&apprshort_atgamma[j],&vtmp1);
      poly_pointwise_montgomery(&tmp2,&apprshort_gamma[j],&vtmp2);
      poly_add(&tmp,&tmp,&tmp2);
      poly_trace65_ntt(&tmp,&tmp);
      poly_shift(&tmp,&tmp,j);
      poly_add(&vprime[i],&vprime[i],&tmp);
    }
  }

  for(i=0;i<R;i++)
    poly_freeze(&vprime[i]);

  return 0;
}
