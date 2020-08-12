#ifndef ROUNDING_H
#define ROUNDING_H

#include <stdint.h>
#include "params.h"

#define power2round_avx NAMESPACE(_power2round_avx)
void power2round_avx(int32_t a1[N], int32_t a0[N], const int32_t a[N]);
#define decompose_avx NAMESPACE(_decompose_avx)
void decompose_avx(int32_t a1[N], int32_t a0[N], const int32_t a[N]);
#define makehint_avx NAMESPACE(_makehint_avx)
void makehint_avx(int32_t h[N], const int32_t a1[N], const int32_t a0[N]);
#define usehint_avx NAMESPACE(_usehint_avx)
void usehint_avx(int32_t b[N], const int32_t a[N], const int32_t h[N]);

#endif
