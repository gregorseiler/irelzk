#include <stdint.h>
#include <string.h>
#include "params.h"
#include "comm.h"
#include "polyvec.h"
#include "poly.h"
#include "consts.h"
#include "product.h"

extern commkey *ck;
extern commrnd y[R];
extern polyvecl yshat[R];
extern polyvecm ymhat[R];

static void autobase(poly r[4], const poly a[4]) {
  poly b[4];

  b[0] = a[0];
  poly_pointwise_montgomery(&b[2],&a[2],&nttx2);
  poly_pointwise_montgomery(&b[1],&a[1],&nttx);
  poly_pointwise_montgomery(&b[3],&a[3],&nttx3);

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
  poly a[4][R], mprime, tmp;
  poly alpha[4], beta[R];
  aes256ctr_ctx aesctx;

  aes256ctr_init(&aesctx,chash,0);
  for(i=0;i<4;i++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&alpha[i],&aesctx);
  }
  for(j=0;j<R;j++) {
    aes256ctr_select(&aesctx,nonce++);
    poly_uniform_preinit(&beta[j],&aesctx);
  }

  /* masking polynomials: a_ij = \sigma^-j(<b_i,y_j>) */
  for(i=0;i<4;i++) {
    for(j=0;j<R;j++) {
      polyvecl_pointwise_acc_montgomery(&a[i][j],&ck->bm[i],&yshat[j]);
      poly_scale_montgomery(&a[i][j],&a[i][j],MONTSQ);
      poly_add(&a[i][j],&a[i][j],&ymhat[j].vec[i]);
    }
    autobase(a[i],a[i]);
  }

  /* g0 */
  memset(v,0,sizeof(poly));
  for(i=0;i<4;i++) {
    for(j=0;j<R;j++) {
      poly_pointwise_montgomery(&tmp,&a[i][j],&a[i][j]);
      poly_pointwise_montgomery(&tmp,&tmp,&alpha[i]);
      poly_pointwise_montgomery(&tmp,&tmp,&beta[j]);
      poly_add(v,v,&tmp);
    }
  }
  poly_scale_montgomery(v,v,MONTSQ);

  for(j=0;j<R;j++) {
    polyvecl_pointwise_acc_montgomery(&tmp,&ck->bm[M-2],&yshat[j]);
    poly_scale_montgomery(&tmp,&tmp,MONTSQ);
    poly_add(&tmp,&tmp,&ymhat[j].vec[M-2]);
    if(j == 1) poly_pointwise_montgomery(&tmp,&tmp,&nttx);
    else if(j == 2) poly_pointwise_montgomery(&tmp,&tmp,&nttx2);
    else if(j == 3) poly_pointwise_montgomery(&tmp,&tmp,&nttx3);
    poly_add(v,v,&tmp);
  }
  poly_freeze(v);

  /* g1 */
  memset(&msg->vec[M-2],0,sizeof(poly));
  for(i=0;i<4;i++) {
    for(j=0;j<N;j++)
      mprime.coeffs[j] = 1 - (msg->vec[i].coeffs[j] << 1);
    for(j=0;j<R;j++) {
      if(j == 0) tmp = mprime;
      else if(j == 1) poly_sigma193_ntt(&tmp,&mprime);
      else if(j == 2) poly_sigma129_ntt(&tmp,&mprime);
      else if(j == 3) poly_sigma65_ntt(&tmp,&mprime);

      for(int k=0;k<N;k++)
        tmp.coeffs[k] *= a[i][j].coeffs[k];

      poly_pointwise_montgomery(&tmp,&tmp,&alpha[i]);
      poly_pointwise_montgomery(&tmp,&tmp,&beta[j]);
      poly_add(&msg->vec[M-2],&msg->vec[M-2],&tmp);
    }
  }
}

int product_verify(poly *v, const uint8_t chash[SYMBYTES], const poly c[R], const commrnd z[R],
                   const comm *tp, const commkey *ckp)
{
  int i,j;
  uint64_t nonce = 0;
  poly tmp, f[R];
  poly chat[R], cfull;
  poly alpha[4], beta[R];
  polyvecl zshat[R];
  polyvecm zmhat[R];
  aes256ctr_ctx aesctx;

  aes256ctr_init(&aesctx,chash,0);
  for(i=0;i<4;i++) {
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
  for(i=0;i<4;i++) {
    for(j=0;j<R;j++) {
      polyvecl_pointwise_acc_montgomery(&f[j],&ckp->bm[i],&zshat[j]);
      poly_pointwise_montgomery(&tmp,&chat[j],&tp->tm.vec[i]);
      poly_sub(&f[j],&f[j],&tmp);
      poly_scale_montgomery(&f[j],&f[j],MONTSQ);
      poly_add(&f[j],&f[j],&zmhat[j].vec[i]);
    }
    autobase(f,f);

    for(j=0;j<R;j++) {
      poly_add(&tmp,&f[j],&cfull);
      poly_pointwise_montgomery(&f[j],&f[j],&tmp);
      poly_pointwise_montgomery(&f[j],&f[j],&alpha[i]);
      poly_pointwise_montgomery(&f[j],&f[j],&beta[j]);
      poly_add(v,v,&f[j]);
    }
  }
  poly_scale_montgomery(v,v,MONTSQ);

  for(j=0;j<R;j++) {
    polyvecl_pointwise_acc_montgomery(&tmp,&ckp->bm[M-2],&zshat[j]);
    poly_scale_montgomery(&tmp,&tmp,MONTSQ);
    poly_add(&tmp,&tmp,&zmhat[j].vec[M-2]);
    if(j == 1) poly_pointwise_montgomery(&tmp,&tmp,&nttx);
    else if(j == 2) poly_pointwise_montgomery(&tmp,&tmp,&nttx2);
    else if(j == 3) poly_pointwise_montgomery(&tmp,&tmp,&nttx3);
    poly_add(v,v,&tmp);
  }
  poly_pointwise_montgomery(&tmp,&cfull,&tp->tm.vec[M-2]);
  poly_scale_montgomery(&tmp,&tmp,MONTSQ);
  poly_sub(v,v,&tmp);
  poly_freeze(v);

  return 0;
}
