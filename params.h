#ifndef PARAMS_H
#define PARAMS_H

//#define ADDITION_PROOF
//#define MULTIPLICATION_PROOF
//#define MULTIPLICATION_PROOF_2
//#define MAIN_PROOF

#define NAMESPACE(s) s

#if defined(ADDITION_PROOF)
#define M 6
#elif defined(MULTIPLICATION_PROOF)
#define M 10
#elif defined(MULTIPLICATION_PROOF_2)
#define M 15
#elif defined(MAIN_PROOF)
#define M 19
#endif

#define N 128
#define Q 1073479681
#define GAMMA1 (1 << 18)
#define GAMMA2 (Q-1)/(1 << 13)
#define D 14
#define BETA 32
#define R 4
#define K 10
#define L 10

#define SYMBYTES 32

#endif
