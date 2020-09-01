#ifndef LINEAR_H
#define LINEAR_H

#include <stdint.h>
#include "params.h"
#include "comm.h"
#include "poly.h"
#include "polyvec.h"

#if defined(ADDITION_PROOF)
void linear(poly vprime[R], poly *h, const polyvecm *msg, uint8_t chash[SYMBYTES]);
int linear_verify(poly vprime[R], const uint8_t chash[SYMBYTES], const poly *h,
                  const poly c[R], const commrnd z[R], const comm *tp, const commkey *ckp);
#elif defined(MULTIPLICATION_PROOF) || defined(MULTIPLICATION_PROOF_2)
void linear_challenge_prehash(const uint8_t chash[SYMBYTES]);
void linear(poly vprime[R], poly *h, const polyvecm *msg, const poly *l, uint8_t chash[SYMBYTES]);
int linear_verify(poly vprime[R], const uint8_t chash[SYMBYTES], const poly *h, const poly *l,
                  const poly c[R], const commrnd z[R], const comm *tp, const commkey *ckp);
#elif defined(MAIN_PROOF)
void linear_first(poly *g);
void linear_challenge_prehash(const uint8_t chash[SYMBYTES], poly A[8][8]);
void linear_last(poly *h, const polyvecm *msg, const poly u[8]);
void linear(poly vprime[R]);
int linear_verify_first(const poly *h);
int linear_verify(poly vprime[R], const poly *h, const poly u[8], const poly c[R], const commrnd z[R],
                  const comm *tp, const commkey *ckp);
#endif

#endif
