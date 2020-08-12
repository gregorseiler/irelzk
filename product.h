#ifndef PRODUCT_H
#define PRODUCT_H

#include <stdint.h>
#include "params.h"
#include "comm.h"
#include "polyvec.h"
#include "poly.h"

void product(poly *v, polyvecm *msg, const uint8_t chash[SYMBYTES]);
int product_verify(poly *v, const uint8_t chash[SYMBYTES], const poly c[4], const commrnd z[4],
                   const comm *tp, const commkey *ckp);

#endif
