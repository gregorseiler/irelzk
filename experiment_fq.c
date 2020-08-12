#include <stdio.h>
#include <stdint.h>

#define Q 1073479681
#define V 262143

int main(void) {
  int32_t i, t, min, max;

  min = (1 << 30) - 1 + (1 << 30);
  max = (1 << 31);
  i = (1 << 31);
  do {
    t  = i >> 30;
    t *= V;
    t  = (i&((1<<30)-1)) + t;

    if(((i%Q)-t)%Q) printf("ERROR\n");
    if(t < min) min = t;
    if(t > max) max = t;
  } while(++i < (1 << 30) - 1 + (1 << 30));
  printf("%d %d\n", min, max);

  return 0;
}
