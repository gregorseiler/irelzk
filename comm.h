#ifndef COMM_H
#define COMM_H

#include "params.h"
#include "polyvec.h"

typedef struct {
  polyveck t0;
  polyvecm tm;
} comm;

typedef struct {
  polyvecl s;
  polyveck e;
  polyvecm em;
} commrnd;

typedef struct {
  polyvecl b0[K];
  polyvecm bt[K];
  polyvecl bm[M];
} commkey;

void expand_commkey(commkey *ck, const uint8_t rho[SYMBYTES]);
void commit(comm *t, commrnd *r, const polyvecm *msg, const commkey *ck);

#endif
