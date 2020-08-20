#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define ntt_avx NAMESPACE(_ntt_avx)
void ntt_avx(int32_t r[N], const int32_t *qdata);
#define invntt_avx NAMESPACE(_invntt_avx)
void invntt_avx(int32_t r[N], const int32_t *qdata);

#endif
