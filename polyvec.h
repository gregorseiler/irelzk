#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

#define POLYVEC_DECLS(POLYVEC_TYPE,POLYVEC_LENGTH) \
  typedef struct { \
    poly vec[POLYVEC_LENGTH]; \
  } POLYVEC_TYPE; \
  void POLYVEC_TYPE##_freeze(POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_add(POLYVEC_TYPE *w, const POLYVEC_TYPE *u, const POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_sub(POLYVEC_TYPE *w, const POLYVEC_TYPE *u, const POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_ntt(POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_invntt(POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_invntt_tomont(POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_pointwise_acc_montgomery(poly *r, const POLYVEC_TYPE *u, const POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_scale_montgomery(POLYVEC_TYPE *v, const POLYVEC_TYPE *u, int32_t s); \
  int POLYVEC_TYPE##_chknorm(const POLYVEC_TYPE *v, uint32_t b); \
  void POLYVEC_TYPE##_power2round(POLYVEC_TYPE *v1, POLYVEC_TYPE *v0, POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_decompose(POLYVEC_TYPE *v1, POLYVEC_TYPE *v0, POLYVEC_TYPE *v); \
  void POLYVEC_TYPE##_makehint(POLYVEC_TYPE *h, const POLYVEC_TYPE *v1, POLYVEC_TYPE *v0); \
  void POLYVEC_TYPE##_usehint(POLYVEC_TYPE *v1, POLYVEC_TYPE *v, const POLYVEC_TYPE *h);

POLYVEC_DECLS(polyveck,K)
POLYVEC_DECLS(polyvecl,L)
POLYVEC_DECLS(polyvecm,M)

#endif
