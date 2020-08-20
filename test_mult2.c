#include <stdint.h>
#include <stdio.h>
#include "cpucycles.h"
#include "speed_print.h"
#include "randombytes.h"
#include "mult.h"
#include "comm.h"

#define NTESTS 500

int main(void) {
  int i, ret;
  uint64_t tsc[NTESTS];
  uint8_t rho[SYMBYTES];
  uint64_t a[2] = {0}, b[2] = {0};
  proof p;
  comm t;
  commrnd r;

  randombytes(rho,SYMBYTES);
  randombytes((uint8_t *)a,16);
  randombytes((uint8_t *)b,16);

  multiplication_proof(&p,&t,&r,rho,a,b);
  ret = multiplication_proof_verify(&p,&t,rho);
  if(ret) printf("FAILURE!\n");

  for(i=0;i<NTESTS;i++) {
    tsc[i] = cpucycles();
    multiplication_proof(&p,&t,&r,rho,a,b);
  }
  print_results("Multiplication Proving: ",tsc,NTESTS);

  for(i=0;i<NTESTS;i++) {
    tsc[i] = cpucycles();
    multiplication_proof_verify(&p,&t,rho);
  }
  print_results("Multiplication Verifying: ",tsc,NTESTS);

  return 0;
}
