#include <stdint.h>
#include <stdio.h>
#include "cpucycles.h"
#include "speed_print.h"
#include "randombytes.h"
#include "addition.h"
#include "comm.h"

#define NTESTS 500

int main(void) {
  int i, ret;
  uint64_t tsc[NTESTS];
  uint8_t rho[SYMBYTES];
  uint64_t a[2], b[2];
  proof p;
  comm t;
  commrnd r;

  randombytes(rho,SYMBYTES);
  randombytes((uint8_t *)a,16);
  randombytes((uint8_t *)b,16);

  addition_proof(&p,&t,&r,rho,a,b);
  ret = addition_proof_verify(&p,&t,rho);
  if(ret) printf("FAILURE!\n");

  for(i=0;i<NTESTS;i++) {
    tsc[i] = cpucycles();
    addition_proof(&p,&t,&r,rho,a,b);
  }
  print_results("Addition Proving: ",tsc,NTESTS);

  for(i=0;i<NTESTS;i++) {
    tsc[i] = cpucycles();
    addition_proof_verify(&p,&t,rho);
  }
  print_results("Addition Verifying: ",tsc,NTESTS);

  return 0;
}
