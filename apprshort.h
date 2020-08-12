#ifndef APPRSHORT_H
#define APPRSHORT_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

void apprshort_challenge_prehash(const uint8_t chash[SYMBYTES]);
void apprshort_challenge_mul(poly *r, const poly *a);
void apprshort_first(poly *yf, const uint8_t seed[SYMBYTES], uint16_t nonce);
int apprshort_last(poly *l, const polyvecm *msg);
int apprshort_verify(const poly *l);

#endif
