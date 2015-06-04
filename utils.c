#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

long
compare(double * t1, double * t2, int n)
{
  long diffs = 0L;
  for(int i = 0; i < n; i++){
    if (fabs(t1[i] - t2[i]) > DOUBLE_COMPARISON_THRESHOLD){
      diffs++;
    }
  }
  return diffs;
}

long
compare_l(long * t1, long * t2, int n)
{
  long diffs = 0L;
  for(int i = 0; i < n; i++){
    if (t1[i] != t2[i]){
      diffs++;
    }
  }
  return diffs;
}

void
init_data_1d(int dimx, double * data)
{
  int i;
  data[0] = 0.0;
  for(i = 1; i < dimx - 1; i ++){
    data[i] = fabs(cos(i) * (1 << 8));
  }
  data[dimx - 1]= 0.0;
}

void
init_data_1d_l(int dimx, long * data)
{
  int i;
  data[0] = 0.0;
  for(i = 1; i < dimx - 1; i ++){
    data[i] = floor(fabs(cos(i) * (1 << 8)));
  }
  data[dimx - 1]= 0.0;
}

void
init_data_2d(int dimx, int dimy, double ** data)
{
  int i,j;
  for(j = 0; j < dimy; j++) {
    data[0][j] = 0.0;
  }
  for(i = 0; i < dimx; i++) {
    data[i][0] = 0.0;
    for(j = 0; j < dimy; j++) {
      data[i][j] = fabs(cos((double) i + j) * (1 << 8));
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
print_test1d_summary(int nruns, int verbose, double total_time,
  struct benchspec bs, double * data_in, double * data_out)
{
  int i;
  if (verbose) {
    printf("\n------------- Input ---------\n");
    for (i = 0; i < DISPLAY_SIZE; i++) {
      printf("%10.3f", data_in[i]);
    }
    printf("Result snapshot: %s\n", KRED);
    for (i = DISPLAY_OFFSET; i < DISPLAY_OFFSET + DISPLAY_SIZE; i++) {
      printf("%10.3f", data_out[i]);
    }
    printf("\n%s----------------------\n", KRESET);
  }
  printf("\n------------- %s ---------\n", bs.name);
  printf("Total time :\t %13f ms\n", (double) total_time * 1000.0);
  printf("Average time :\t %13f ms\n\n",
    (double) (total_time * 1000.0 / (nruns)));

}


void
print_test1d_l_summary(int nruns, int verbose, double total_time,
  struct benchspec1d_l bs, long * data_in, long * data_out)
{
  int i;
  if (verbose) {
    printf("\n------------- Input ---------\n");
    for (i = 0; i < DISPLAY_SIZE; i++) {
      printf(" %4li ", data_in[i]);
    }
    printf("Result snapshot: %s\n", KRED);
    for (i = DISPLAY_OFFSET; i < DISPLAY_OFFSET + DISPLAY_SIZE; i++) {
      printf(" %4li ", data_out[i]);
    }
    printf("\n------------- %s ---------\n", bs.name);
  }
  printf("\n%s----------------------\n", KRESET);
  printf("Total time :\t %13f ms\n", (double) total_time * 1000.0);
  printf("Average time :\t %13f ms\n\n",
    (double) (total_time * 1000.0 / (nruns)));

}

void
print_test2d_summary(int nruns, int verbose, double total_time,
  struct benchspec2d bs, double ** data_in, double ** data_out)
{
  int i,j;
  if (verbose) {
    printf("\n------------- Input ---------\n");
    for (i = DISPLAY_OFFSET; i < DISPLAY_OFFSET + DISPLAY_SIZE; i++) {
      for (j = DISPLAY_OFFSET; j < DISPLAY_OFFSET + DISPLAY_SIZE; j++) {
        printf("%10.3f", data_in[i][j]);
      }
      printf("\n");
    }
    printf("Result snapshot: \n");
    for (i = DISPLAY_OFFSET; i < DISPLAY_OFFSET + DISPLAY_SIZE; i++) {
      for (j = DISPLAY_OFFSET; j < DISPLAY_OFFSET + DISPLAY_SIZE; j++) {
        printf("%10.3f", data_out[i][j]);
      }
      printf("\n");
    }
    printf("\n----------------------\n");
  }
  printf("\n------------- %s ---------\n", bs.name);
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
