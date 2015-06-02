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
init_data_1d(int dimx, double * data)
{
  int i;
  for(i = 0; i < dimx; i ++){
    data[i] = cos(i) * (1 << 8);
  }
}

void
init_data_2d(int dimx, int dimy, double ** data)
{
  int i,j;
  for(i = 0; i < dimx; i++) {
    for(j = 0; j < dimy; j++) {
      data[i][j] = cos((double) i + j) * (1 << 8);
    }
  }
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
print_runscores(int nruns, struct benchscore * bsc)
{
  int i;
  for (i = 0; i < nruns; i++){
    printf("%s %s %s : Run %i ...", KBLU, bsc[i].name, KRESET, i + 1);
    printf("\t\t %13f ms\n", bsc[i].wallclock * 1000.0 );
  }
}

void
print_test2d_summary(int nruns, double total_time, struct benchspec2d bs,
  double **data_out)
{
  int i,j;
  printf("\n------------- %s ---------\n", bs.name);
  printf("Result snapshot: \n");
  for (i = 0; i < DISPLAY_SIZE; i++) {
    for (j = 0; j < DISPLAY_SIZE; j++) {
      printf("%10.3f", data_out[i][j]);
    }
  }
  printf("\n----------------------\n");
  printf("Total time :\t %13f ms\n", (double) total_time * 1000.0);
  printf("Average time :\t %13f ms\n\n",
    (double) (total_time * 1000.0 / (nruns)));

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
