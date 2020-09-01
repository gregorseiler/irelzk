#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "params.h"
#include "comm.h"
#include "polyvec.h"
#include "poly.h"
#include "consts.h"
#include "product.h"

extern const polyvecm g[R];

static inline void poly_sigmainv_ntt(poly *r, const poly *a, int i) {
  if(i == 0) *r = *a;
  else if(i == 1) poly_sigma193_ntt(r,a);
  else if(i == 2) poly_sigma129_ntt(r,a);
  else if(i == 3) poly_sigma65_ntt(r,a);
}

static inline void poly_shift(poly *r, const poly *a, int i) {
  if(i == 0) *r = *a;
  else if(i == 1) poly_pointwise_montgomery(r,a,&nttx);
  else if(i == 2) poly_pointwise_montgomery(r,a,&nttx2);
  else if(i == 3) poly_pointwise_montgomery(r,a,&nttx3);
}

static void autobase(poly r[4], const poly *a, size_t off) {
  poly b[4];

  b[0] = a[0*off];
  poly_pointwise_montgomery(&b[1],&a[1*off],&nttx);
  poly_pointwise_montgomery(&b[2],&a[2*off],&nttx2);
  poly_pointwise_montgomery(&b[3],&a[3*off],&nttx3);

  poly_add(&r[0],&b[0],&b[2]);
  poly_add(&r[1],&b[1],&b[3]);
  poly_sub(&r[2],&r[0],&r[1]);
  poly_add(&r[0],&r[0],&r[1]);

  poly_sub(&b[1],&b[1],&b[3]);
  poly_pointwise_montgomery(&b[1],&b[1],&nttx64);
  poly_sub(&b[0],&b[0],&b[2]);
  poly_add(&r[1],&b[0],&b[1]);
  poly_sub(&r[3],&b[0],&b[1]);

  poly_sigma193_ntt(&r[1],&r[1]);
  poly_sigma129_ntt(&r[2],&r[2]);
  poly_sigma65_ntt(&r[3],&r[3]);
}

void product(poly *v, polyvecm *msg, const uint8_t chash[SYMBYTES]) {
  int i,j;
  uint64_t nonce = 0;
  poly a[R], mprime, tmp, tmp2, tmp3;
  poly alpha[16], beta[R];
  aes256ctr_ctx aesctx;

  //FIXME: Only once
  aes256ctr_init(&aesctx,chash,0);
  for(i=0;i<16;i++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&alpha[i],&aesctx);
  }
  for(j=0;j<R;j++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&beta[j],&aesctx);
  }

  memset(v,0,sizeof(poly));
  memset(&msg->vec[M-3],0,sizeof(poly));
  memset(&msg->vec[M-2],0,sizeof(poly));
  for(i=0;i<16;i++) {
    /* a_j = \sigma^-j(<b_i,y_j>) */
    autobase(a,&g[0].vec[i],M);

    /* garbage */
    for(j=0;j<N;j++)
      mprime.coeffs[j] = 3*msg->vec[i].coeffs[j];
    for(j=0;j<R;j++) {
      poly_pointwise_montgomery(&tmp,&a[j],&alpha[i]);
      poly_pointwise_montgomery(&tmp,&tmp,&beta[j]);
      poly_sigmainv_ntt(&tmp2,&mprime,j);

      poly_sigmainv_ntt(&tmp3,&msg->vec[i],j);
      for(int k=0;k<N;k++)
        tmp3.coeffs[k] = tmp2.coeffs[k]*tmp3.coeffs[k] - 1;
      poly_pointwise_montgomery(&tmp3,&tmp,&tmp3);
      poly_add(&msg->vec[M-2],&msg->vec[M-2],&tmp3);

      poly_pointwise_montgomery(&tmp,&tmp,&a[j]);
      poly_pointwise_montgomery(&tmp2,&tmp,&tmp2);
      poly_sub(&msg->vec[M-3],&msg->vec[M-3],&tmp2);

      poly_pointwise_montgomery(&tmp,&tmp,&a[j]);
      poly_add(v,v,&tmp);
    }
  }

  poly_scale_montgomery(v,v,65810308); // mont^3
  poly_scale_montgomery(&msg->vec[M-3],&msg->vec[M-3],65810308);
  poly_scale_montgomery(&msg->vec[M-2],&msg->vec[M-2],MONTSQ);
  for(i=0;i<R;i++) {
    poly_shift(&tmp,&g[i].vec[M-2],i);
    poly_add(&msg->vec[M-3],&msg->vec[M-3],&tmp);
  }
  for(i=0;i<R;i++) {
    poly_shift(&tmp,&g[i].vec[M-3],i);
    poly_add(v,v,&tmp);
  }
  poly_freeze(v);
}

int product_verify(poly *v, const uint8_t chash[SYMBYTES], const poly c[R], const commrnd z[R],
                   const comm *tp, const commkey *ckp)
{
  int i,j;
  uint64_t nonce = 0;
  poly f[R], tmp;
  poly alpha[16], beta[R];
  poly chat[R], cfull;
  polyvecl zshat[R];
  polyvecm zmhat[R];
  aes256ctr_ctx aesctx;

  aes256ctr_init(&aesctx,chash,0);
  for(i=0;i<16;i++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&alpha[i],&aesctx);
  }
  for(j=0;j<R;j++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&beta[j],&aesctx);
  }

  for(j=0;j<R;j++) {
    zshat[j] = z[j].s;
    zmhat[j] = z[j].em;
    chat[j] = c[j];
    polyvecl_ntt(&zshat[j]);
    polyvecm_ntt(&zmhat[j]);
    poly_ntt(&chat[j]);
  }

  poly_pointwise_montgomery(&tmp,&chat[1],&nttx);
  poly_add(&cfull,&chat[0],&tmp);
  poly_pointwise_montgomery(&tmp,&chat[2],&nttx2);
  poly_add(&cfull,&cfull,&tmp);
  poly_pointwise_montgomery(&tmp,&chat[3],&nttx3);
  poly_add(&cfull,&cfull,&tmp);

  memset(v,0,sizeof(poly));
  for(i=0;i<16;i++) {
    for(j=0;j<R;j++) {
      polyvecl_pointwise_acc_montgomery(&f[j],&ckp->bm[i],&zshat[j]);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[i]);
      poly_sub(&f[j],&f[j],&tmp);
      poly_scale_montgomery(&f[j],&f[j],MONTSQ);
      poly_add(&f[j],&f[j],&zmhat[j].vec[i]);
    }
    autobase(f,f,1);

    for(j=0;j<R;j++) {
      poly_add(&tmp,&f[j],&cfull);
      poly_pointwise_montgomery(&tmp,&tmp,&f[j]);
      poly_sub(&f[j],&f[j],&cfull);
      poly_pointwise_montgomery(&tmp,&tmp,&f[j]);
      poly_pointwise_montgomery(&tmp,&tmp,&alpha[i]);
      poly_pointwise_montgomery(&tmp,&tmp,&beta[j]);
      poly_add(v,v,&tmp);
    }
  }

  poly_scale_montgomery(v,v,MONTSQ);
  for(j=0;j<R;j++) {
    polyvecl_pointwise_acc_montgomery(&tmp,&ckp->bm[M-2],&zshat[j]);
    poly_scale_montgomery(&tmp,&tmp,MONTSQ);
    poly_add(&tmp,&tmp,&zmhat[j].vec[M-2]);
    poly_shift(&tmp,&tmp,j);
    poly_pointwise_montgomery(&tmp,&tmp,&cfull);
    poly_add(v,v,&tmp);
  }
  poly_pointwise_montgomery(&tmp,&cfull,&tp->tm.vec[M-2]);
  poly_pointwise_montgomery(&tmp,&cfull,&tmp);
  poly_scale_montgomery(&tmp,&tmp,MONTSQ);
  poly_sub(v,v,&tmp);

  poly_scale_montgomery(v,v,MONTSQ);
  for(j=0;j<R;j++) {
    polyvecl_pointwise_acc_montgomery(&tmp,&ckp->bm[M-3],&zshat[j]);
    poly_scale_montgomery(&tmp,&tmp,MONTSQ);
    poly_add(&tmp,&tmp,&zmhat[j].vec[M-3]);
    poly_shift(&tmp,&tmp,j);
    poly_add(v,v,&tmp);
  }
  poly_pointwise_montgomery(&tmp,&cfull,&tp->tm.vec[M-3]);
  poly_scale_montgomery(&tmp,&tmp,MONTSQ);
  poly_sub(v,v,&tmp);

  poly_freeze(v);

  return 0;
}
