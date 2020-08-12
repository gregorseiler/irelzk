#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "aes256ctr.h"
#include "params.h"

typedef struct {
  int32_t coeffs[N];
} poly __attribute__((aligned(32)));

#define poly_uniform_preinit NAMESPACE(_poly_uniform_preinit)
void poly_uniform_preinit(poly *r, aes256ctr_ctx *state);
#define poly_uniform NAMESPACE(_poly_uniform)
void poly_uniform(poly *r, const uint8_t seed[SYMBYTES], uint16_t nonce);
#define poly_trinary_preinit NAMESPACE(_poly_trinary_preinit)
void poly_trinary_preinit(poly *r, aes256ctr_ctx *state);
#define poly_trinary NAMESPACE(_poly_trinary)
void poly_trinary(poly *r, const uint8_t seed[SYMBYTES], uint16_t nonce);
#define poly_uniform_gamma_preinit NAMESPACE(_poly_uniform_gamma_preinit)
void poly_uniform_gamma_preinit(poly *r, aes256ctr_ctx *state);
#define poly_uniform_gamma NAMESPACE(_poly_uniform_gamma)
void poly_uniform_gamma(poly *r, const uint8_t seed[SYMBYTES], uint16_t nonce);
#define poly_reduce NAMESPACE(_poly_reduce)
void poly_reduce(poly *r);
#define poly_freeze NAMESPACE(_poly_freeze)
void poly_freeze(poly *r);
#define poly_csubq NAMESPACE(_poly_csubq)
void poly_csubq(poly *r);
#define poly_caddq NAMESPACE(_poly_caddq)
void poly_caddq(poly *r);
#define poly_neg NAMESPACE(_poly_neg)
void poly_neg(poly *r, const poly *a);
#define poly_add NAMESPACE(_poly_add)
void poly_add(poly *r, const poly *a, const poly *b);
#define poly_sub NAMESPACE(_poly_sub)
void poly_sub(poly *r, const poly *a, const poly *b);
#define poly_ntt NAMESPACE(_poly_ntt)
void poly_ntt(poly *r);
#define poly_invntt NAMESPACE(_poly_invntt)
void poly_invntt(poly *r);
#define poly_invntt_tomont NAMESPACE(_poly_invntt_tomont)
void poly_invntt_tomont(poly *r);
#define poly_pointwise_montgomery NAMESPACE(_poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *r, const poly *a, const poly *b);
#define poly_scale_montgomery NAMESPACE(_poly_scale_montgomery)
void poly_scale_montgomery(poly *r, const poly *a, int32_t s);
#define poly_ntt2 NAMESPACE(_poly_ntt2)
void poly_ntt2(poly *r);
#define poly_invntt2 NAMESPACE(_poly_invntt2)
void poly_invntt2(poly *r);
#define poly_chknorm NAMESPACE(_poly_chknorm)
int poly_chknorm(const poly *a, uint32_t b);
#define poly_rotate NAMESPACE(_poly_rotate)
void poly_rotate(poly *r, const poly *a, int k);
#define poly_sigma NAMESPACE(_poly_sigma)
void poly_sigma(poly *r, const poly *a, int k);
#define poly_sigma65_ntt NAMESPACE(_poly_sigma65_ntt)
void poly_sigma65_ntt(poly *r, const poly *a);
#define poly_sigma129_ntt NAMESPACE(_poly_sigma129_ntt)
void poly_sigma129_ntt(poly *r, const poly *a);
#define poly_sigma193_ntt NAMESPACE(_poly_sigma193_ntt)
void poly_sigma193_ntt(poly *r, const poly *a);
#define poly_trace65_ntt NAMESPACE(_poly_trace65_ntt)
void poly_trace65_ntt(poly *r, const poly *a);
#define poly_power2round NAMESPACE(_poly_power2round)
void poly_power2round(poly *a1, poly *a0, poly *a);
#define poly_decompose NAMESPACE(_poly_decompose)
void poly_decompose(poly *a1, poly *a0, poly *a);
#define poly_makehint NAMESPACE(_poly_makehint)
void poly_makehint(poly *h, const poly *a1, poly *a0);
#define poly_usehint NAMESPACE(_poly_usehint)
void poly_usehint(poly *b1, poly *a, const poly *h);

#endif
