#ifndef MULT_H
#define MULT_H

#include <stdint.h>
#include "params.h"
#include "comm.h"

typedef struct {
  poly h;
  poly c[4];
  commrnd z[4];
} proof;

void main_proof(proof *p, comm *t, commrnd *r, const uint8_t rho[SYMBYTES],
                poly At[8][8], const poly u[8], const poly s[16]);
int main_proof_verify(const proof *p, const comm *t, const uint8_t rho[SYMBYTES],
                      poly At[8][8], const poly u[8]);

#endif
