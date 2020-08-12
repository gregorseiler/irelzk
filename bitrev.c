#include <stdio.h>
#include <stdint.h>

static uint8_t bitrev7(uint8_t a) {
  uint8_t t;

  t  = (a & 1) << 6;
  t |= (a & 2) << 4;
  t |= (a & 4) << 2;
  t |= (a & 8) << 0;
  t |= (a & 16) >> 2;
  t |= (a & 32) >> 4;
  t |= (a & 64) >> 6;
  return t;
}

int main(void) {
  unsigned int i, j, k;
  uint8_t t;

  printf("%hu, ", bitrev7(1));
  printf("%hu, ", bitrev7(2));
  printf("%hu, ", bitrev7(3));

  for(i = 0; i < 4; ++i) {
    for(j = 4; j <= 32; j <<= 1)
      for(k = 0; k < j/4; ++k)
        printf("%hu, ", bitrev7(j + (j/4)*i + k));

    for(k = 0; k < 8; ++k)
      printf("%hu, ", bitrev7(64 + 16*i + 2*k));
    for(k = 0; k < 8; ++k)
      printf("%hu, ", bitrev7(64 + 16*i + 2*k+1));
  }
  printf("\n");
}
