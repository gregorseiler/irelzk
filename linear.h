#ifndef LINEAR_H
#define LINEAR_H

#include <stdint.h>
#include "params.h"
#include "comm.h"
#include "poly.h"
#include "polyvec.h"

void linear_challenge_prehash(const uint8_t chash[SYMBYTES]);
#if defined(ADDITION_PROOF)
void linear(poly *h, poly vprime[R], const polyvecm *msg, uint8_t chash[SYMBYTES]);
int linear_verify(poly vprime[R], const uint8_t chash[SYMBYTES], const poly *h,
                  const poly c[R], const commrnd z[R], const comm *tp, const commkey *ckp);
#elif defined(MULTIPLICATION_PROOF) || defined(MULTIPLICATION_PROOF_2)
void linear(poly *h, poly vprime[R], const polyvecm *msg, const poly *l, uint8_t chash[SYMBYTES]);
int linear_verify(poly vprime[R], const uint8_t chash[SYMBYTES], const poly *h, const poly *l,
                  const poly c[R], const commrnd z[R], const comm *tp, const commkey *ckp);
#endif

#endif
