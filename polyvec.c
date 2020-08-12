#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"

typedef POLYVEC_TYPE polyvec;
#define _CAT(a,b) a##b
#define CAT(a,b) _CAT(a,b)
#define POLYVEC_NAMESPACE(s) CAT(POLYVEC_TYPE,s)
#define polyvec_freeze POLYVEC_NAMESPACE(_freeze)
#define polyvec_add POLYVEC_NAMESPACE(_add)
#define polyvec_sub POLYVEC_NAMESPACE(_sub)
#define polyvec_ntt POLYVEC_NAMESPACE(_ntt)
#define polyvec_invntt POLYVEC_NAMESPACE(_invntt)
#define polyvec_invntt_tomont POLYVEC_NAMESPACE(_invntt_tomont)
#define polyvec_pointwise_acc_montgomery POLYVEC_NAMESPACE(_pointwise_acc_montgomery)
#define polyvec_scale_montgomery POLYVEC_NAMESPACE(_scale_montgomery)
#define polyvec_chknorm POLYVEC_NAMESPACE(_chknorm)
#define polyvec_power2round POLYVEC_NAMESPACE(_power2round)
#define polyvec_decompose POLYVEC_NAMESPACE(_decompose)
#define polyvec_makehint POLYVEC_NAMESPACE(_makehint)
#define polyvec_usehint POLYVEC_NAMESPACE(_usehint)

void polyvec_freeze(polyvec *v) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_freeze(&v->vec[i]);
}

void polyvec_add(polyvec *w, const polyvec *u, const polyvec *v) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyvec_sub(polyvec *w, const polyvec *u, const polyvec *v) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyvec_ntt(polyvec *v) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_ntt(&v->vec[i]);
}

void polyvec_invntt(polyvec *v) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_invntt(&v->vec[i]);
}

void polyvec_invntt_tomont(polyvec *v) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_invntt_tomont(&v->vec[i]);
}

void polyvec_pointwise_acc_montgomery(poly *r, const polyvec *u, const polyvec *v) {
  int i;
  poly t;

  poly_pointwise_montgomery(r, &u->vec[0], &v->vec[0]);
  for(i=1;i<POLYVEC_LENGTH;i++) {
    poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
    poly_add(r, r, &t);
  }
}

void polyvec_scale_montgomery(polyvec *v, const polyvec *u, int32_t s) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_scale_montgomery(&v->vec[i], &u->vec[i], s);
}

int polyvec_chknorm(const polyvec *v, uint32_t b) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    if(poly_chknorm(&v->vec[i], b))
      return 1;

  return 0;
}

void polyvec_power2round(polyvec *v1, polyvec *v0, polyvec *v) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_power2round(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

void polyvec_decompose(polyvec *v1, polyvec *v0, polyvec *v) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_decompose(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

void polyvec_makehint(polyvec *h, const polyvec *v1, polyvec *v0) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_makehint(&h->vec[i], &v1->vec[i], &v0->vec[i]);
}

void polyvec_usehint(polyvec *v1, polyvec *v, const polyvec *h) {
  int i;
  for(i=0;i<POLYVEC_LENGTH;i++)
    poly_usehint(&v1->vec[i], &v->vec[i], &h->vec[i]);
}

