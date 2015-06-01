#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

int
compare(double * t1, double * t2, int n)
{
  for(int i = 0; i < n; i++){
    if(fabs(t1[i] - t2[i]) > DOUBLE_COMPARISON_THRESHOLD){
      return 0;
    }
  }
  return 1;
}

void
print_benchspecs(int n, struct benchspec * benchmarks)
{
  int i;
  for (i = 0; i < n; i++) {
    printf("Bench:%s %s %s\n Iterations :\t%s%i%s \tProblem size :\t%s%li%s\n",
      KRED, benchmarks[i].name, KRESET,
      KCYN, benchmarks[i].iters, KRESET,
      KCYN, benchmarks[i].size , KRESET);
  }
}


void
swap(void *a, void *b, size_t size)
{
  enum { threshold = (1 << 7) };
  if (size <= threshold) {
    char temp[threshold];

    memcpy(temp, b,    size);
    memcpy(b,    a,    size);
    memcpy(a,    temp, size);
  }
  else {
    void* temp = aligned_alloc(CACHE_LINE_SIZE, size);

    memcpy(temp, b,    size);
    memcpy(b,    a,    size);
    memcpy(a,    temp, size);

    free(temp);
  }
}
