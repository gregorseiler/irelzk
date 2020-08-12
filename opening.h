#ifndef OPENING_H
#define OPENING_H

#include <stdint.h>
#include "params.h"
#include "comm.h"
#include "polyvec.h"
#include "poly.h"

void opening_init(const polyveck *t0lowp, const commrnd *rp, const commkey *ckp);
void opening_first(polyveck w1[4], const uint8_t seed[SYMBYTES], uint16_t nonce);

void challenge_prehash(poly c[4], const uint8_t chash[N/4]);
void challenge(poly c[4], const polyveck w[4]);

int opening_last(commrnd z[4], const poly c[4], const polyveck w1[R]);

int opening_verify_first(polyveck w1[4], const poly c[4], const commrnd z[4],
                         const comm *tp, const commkey *ckp);
int opening_verify_last(const poly c[4], const uint8_t chash[N/4]);

#endif
