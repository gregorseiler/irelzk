#include <stdio.h>
#include <stdint.h>

int main(void) {
  int32_t i,a0,a1;
  const int32_t q = 1073479681;
  const int32_t alpha = (q-1)/(1 << 12);

  for(i=-(q-1)/2;i<=(q-1)/2;i++) {
    a1 = (i + (i >> 12) + (i >> 24) + (1 << 17)) >> 18;
    a0 = i - a1*alpha;
    if(a0 > alpha/2)
      printf("ERROR4: %d\n", i);

    if(a0 <= -alpha/2)
      printf("ERROR5: %d\n", i);
  }

  return 0;
}
