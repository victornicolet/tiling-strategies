#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

long
compare(int n, double *t1, double *t2)
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
compare_fast(int n, double *t1, double *t2)
{
  long diffs = 0L;
  #pragma omp parallel for schedule(static) reduction(+ : diffs)
  for(int i = 0; i < n; i++){
    if (fabs(t1[i] - t2[i]) > DOUBLE_COMPARISON_THRESHOLD){
      diffs++;
    }
  }
  return diffs;
}

void print_check(int diffs, int pb_size, double *out_mx, double *ref_out)
{
  int i;
  int *diff_rec;

  diff_rec = malloc(pb_size * sizeof(*diff_rec));

  find_diffs(pb_size, out_mx, ref_out, diff_rec);
  (void)printf("Output not correct : %i/%i -(%4.3f)."
    "\n",
    diffs,
    pb_size,
    (double) diffs / pb_size);
  (void)printf("Output from this version\n");

  for (i = 0; i < DISPLAY_SIZE; i++) {
    (void)printf("%10.3f", out_mx[i]);
  }
  (void)printf("\nOutput from reference version :\n");

  for (i = 0; i < DISPLAY_SIZE; i++) {
    (void)printf("%10.3f", ref_out[i]);
  }

  (void)printf("\n");
  (void)printf("\nDiff record snapshot :\n");

  for (i = 0; i < min(DISPLAY_SIZE, diffs); i++) {
    (void)printf("%10i", diff_rec[i]);
  }

  (void)printf("\n");

  free(diff_rec);
}

void find_diffs(int n, double *a, double *b, int * diff_rec){
        /* Wno-unused */
  int i,j;

  j = 0;
  for (i = 0; i < n; i++) {
    if (fabs(a[i]-b[i]) > DOUBLE_COMPARISON_THRESHOLD) {
      diff_rec[j++] = i;
    }
  }
}

long
compare_l(long *t1, long *t2, int n)
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
init_data_1d(int dimx, double *data)
{
  int i;
  data[0] = 0.0;
  for(i = 1; i < dimx - 1; i ++){
    data[i] = fabs(cos(i) * (1 << 8));
  }
  data[dimx - 1]= 0.0;
}

void
init_data_1d_l(int dimx, long *data)
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
print_test1d_summary(int nruns, int verbose, double total_time,
  struct benchspec bs, double * data_in, double * data_out)
{
  int i;
  if (verbose) {
    printf("\n------------- Input ---------\n");
    for (i = DISPLAY_OFFSET; i < DISPLAY_OFFSET + DISPLAY_SIZE; i++) {
      printf("%10.3f", data_in[i]);
    }
    printf("\nResult snapshot: %s\n", KRED);
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
    printf("\nResult snapshot: %s\n", KRED);
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
    printf("\nResult snapshot: \n");
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
swap(double **a, double **b)
{
    double *tmp = *a;
    *a = *b;
    *b = tmp;
}

void
why_fopen(int err_no)
{
  switch(errno){
    case EACCES:
      fprintf(stderr, "You don't have access to this files\n");
      break;
    case EFAULT:
      fprintf(stderr, "Pathname points outside your accessible address space\n");
      break;
    case ENOENT:
      fprintf(stderr, "Pathname refers to a non-existent directory !\n");
    case ENOMEM:
      fprintf(stderr, "Insufficient kernel memory available.\n");
      break;
    default:
      fprintf(stderr, "Errno not yet impl. here\n");
  }
}
