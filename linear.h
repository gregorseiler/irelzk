#ifndef LINEAR_H
#define LINEAR_H

#include <stdint.h>
#include "params.h"
#include "comm.h"
#include "poly.h"
#include "polyvec.h"

void linear_challenge_prehash(const uint8_t chash[SYMBYTES]);
#if defined(ADDITION_PROOF)
void linear(poly *h, poly vprime[4], const polyvecm *msg, uint8_t chash[SYMBYTES]);
int linear_verify(poly vprime[4], const uint8_t chash[SYMBYTES], const poly *h,
                  const poly c[4], const commrnd z[4], const comm *tp, const commkey *ckp);
#elif defined(MULTIPLICATION_PROOF)
void linear(poly *h, poly vprime[4], const polyvecm *msg, const poly *l, uint8_t chash[SYMBYTES]);
int linear_verify(poly vprime[4], const uint8_t chash[SYMBYTES], const poly *h, const poly *l,
                  const poly c[4], const commrnd z[4], const comm *tp, const commkey *ckp);
#endif

#endif
