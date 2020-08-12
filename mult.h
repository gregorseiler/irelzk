#ifndef MULT_H
#define MULT_H

#include <stdint.h>
#include "params.h"
#include "comm.h"

typedef struct {
  poly h;
  poly l;
  poly c[4];
  commrnd z[4];
} proof;

void multiplication_proof(proof *p, comm *t, commrnd *r, const uint8_t rho[SYMBYTES],
                          const uint64_t a[2], const uint64_t b[2]);
int multiplication_proof_verify(const proof *p, const comm *t, const uint8_t rho[SYMBYTES]);

#endif
