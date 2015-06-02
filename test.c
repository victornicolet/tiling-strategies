#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>

#include "utils.h"

#include "benchmarks/jacobi1d.h"
#include "benchmarks/jacobi2d.h"

/* Problem size (in space) between 2 ^ MIN_POW and 2 ^ MAX_POW */
#define MIN_POW 10
#define MAX_POW 16


//static FILE * csv_file;

double test1d(int, int, int, struct benchspec);
double test2d(int, int, int, int, struct benchspec2d);
struct args_dimt get2dargs(int, int, int, struct benchspec2d);
struct args_dimt get1dargs(int, int, struct benchspec);
void usage(int, int, char **, struct benchspec *, struct benchspec2d *);

int main(int argc, char ** argv){

  struct benchspec benchmarks[] = {
    {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap, check_default,
      1 >> 15, 1 >> 5, 1},
    {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive, check_default,
      1 >> 15, 1 >> 5, 1},
    {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test, check_default,
      1 >> 15, 1 >> 5, 1},
    {"JACOBI1D_SWAP_SEQ", djbi1d_sequential, check_default,
      1 >> 15, 1 >> 5, 1},
    {"JACOBI1D_HALF_DIAMONDS", djbi1d_half_diamonds_test, check_low_iter,
      1 >> 15, 1 >> 5, 1},
  };

  struct  benchspec2d benchmarks2d [] = {
    {"JACOBI2D_HALF_DIAMONDS", djbi2d_half_diamonds, check2d_default,
      1 >> 10, 1 >> 10, 1 >> 5},
    {"JACOBI2D_SEQ ", djbi2d_seq, check2d_default,
      1 >> 10, 1 >> 10, 1 >> 5},

  };

  int nbench = sizeof(benchmarks) / sizeof(struct benchspec);
  int nbench2d = sizeof(benchmarks2d) / sizeof(struct benchspec2d);

  /*---------- Parameters section -----------*/

  if(argc < 2){
    usage(nbench, nbench2d, argv, benchmarks, benchmarks2d);
    return  0;
  }

  char * benchmask = argv[1];
  int nruns = 0;
  int maskl;

  if(argc == 3){
    nruns = atoi(argv[2]);
  }
  if((maskl = strlen(benchmask)) > nbench + nbench2d){
    printf("Error : not a valid mask ! Your mask must be %i bits long\n",
      nbench + nbench2d);
    return -1;
  }

  int dimx = 0, dimy = 0, dimt = 0;

  #ifdef DEBUG
    dimt = DEBUG_ITERS ;
    dimx = DEBUG_SIZE;
    dimy = DEBUG_SIZE;
  #else
    if(argc >= 6){
      dimt = atoi(argv[3]);
      dimx = atoi(argv[4]);
      dimy = atoi(argv[5]);
    }
  #endif
  /* -------------------------------------- */

  double exec_time = 0.0, prev_xctime = -1.0;
  for(int bs = 0; bs < maskl; bs++){
    exec_time = 0.0;
    if (benchmask[bs] == '1' && bs < nbench) {
      printf("Execute test for %s ...\n", benchmarks[bs].name);
      exec_time = test1d(nruns, dimx, dimt, benchmarks[bs]);
    } else if (benchmask[bs] == '1'){
      printf("Execute test for %s ...\n", benchmarks2d[bs - nbench].name);
      exec_time = test2d(nruns, dimx, dimy, dimt, benchmarks2d[bs - nbench]);
    }
    if(exec_time < prev_xctime){
      printf("%s %10.3f better !%s\n", KRED, exec_time, KRESET);
    }
  }
  return 0;
}

double
test2d(int nruns, int dimx, int dimy, int dimt, struct benchspec2d benchmark)
{
  int i;
  double t_accu = 0.0;

  double **data_in, **data_out;

  struct args_dimt args = get2dargs(dimx, dimy, dimt, benchmark);
  struct benchscore scores[nruns];

  data_in = alloc_double_mx(args.width, args.height);
  data_out = alloc_double_mx(args.width, args.height);
  init_data_2d(args.width, args.height, data_in);

  for (i = 0; i <= nruns; i ++) {
    printf("%i,", i);
    if (i == 0) {
      benchmark.variant(args, data_in, &scores[0], data_out);
    } else {
      benchmark.variant(args, data_in, &scores[i - 1], data_out);
      scores[i - 1].name = benchmark.name;
      t_accu += scores[i - 1].wallclock;
    }
  }
  printf("Done !\n");
  print_runscores(nruns, scores);
  print_test2d_summary(nruns, t_accu, benchmark, data_out);
  free_mx(data_out, dimx);
  free_mx(data_in, dimx);
  return 0.0;
}

double
test1d(int nruns, int dimx, int dimt, struct benchspec benchmark)
{
  int i;
  double t_accu;
  struct args_dimt args = get1dargs(dimx, dimt, benchmark);

  double * data_in = (double *) aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(*data_in) * args.width);
  double * data_out = (double *) aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(*data_out) * args.width);

  init_data_1d(args.width, data_in);

  struct benchscore scores[nruns];

  t_accu = 0.0;
  for (i = 0; i <= nruns; i ++) {
    printf("%i,", i);
    if (i == 0) {
      benchmark.variant(args, data_in, data_out, &scores[0]);
    } else {
      benchmark.variant(args, data_in, data_out, &scores[i - 1]);
      scores[i - 1].name = benchmark.name;
      t_accu += scores[i - 1].wallclock;
    }
  }
  free(data_in);
  free(data_out);
  return 0.0;
}

struct args_dimt
get2dargs(int dimx, int dimy, int dimt, struct benchspec2d bs)
{
  if(dimx == 0 || dimy == 0 || dimt == 0){
    dimx = bs.width;
    dimy = bs.height;
    dimt = bs.iters;
  }
  struct args_dimt res;
  res.width = dimx;
  res.height = dimy;
  res.iters = dimt;
  return res;
}

struct args_dimt
get1dargs(int dimx, int dimt, struct benchspec bs)
{
  if(dimx == 0 || dimt == 0){
    dimx = bs.size;
    dimt = bs.iters;
  }
  struct args_dimt res;
  res.width = dimx;
  res.height = 0;
  res.iters = dimt;
  return res;
}

void
usage(int nbs, int nbs2d, char ** argv, struct benchspec * bs,
  struct benchspec2d * bs2d)
{
  int i;
  printf("Usage : %s <Mask (%i bit wide) > <Test runs / alg.>\n", argv[0],
    nbs + nbs2d);
  printf("Available benchmarks  :\n");
  printf("%s 1D Benchmarks %s\n", KBLU, KRESET);
  for(i = 0; i < nbs; i++){
    printf("%i - %s\n", i, bs[i].name);
  }
  printf("%s 2D Benchmarks %s\n", KBLU, KRESET);
  for(i = 0; i < nbs2d; i++){
    printf("%i - %s\n", i, bs2d[i].name);
  }
}
