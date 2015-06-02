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
struct args_dimt get2dargs(int, int, int, double **, struct benchspec2d);

int main(int argc, char ** argv){

  int i;

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
      1 >> 15, 1 >> 15, 1 >> 5},
    {"JACOBI2D_SEQ ", djbi2d_seq, check2d_default,
      1 >> 15, 1 >> 15, 1 >> 5},

  };

  int nbench = sizeof(benchmarks) / sizeof(struct benchspec);
  int nbench2d = sizeof(benchmarks2d) /sizeof(struct benchspec2d);

  /*---------- Parameters section -----------*/

  if(argc < 2){
    printf("Usage : %s <Mask (%i bit wide) > <Test runs / alg.>\n", argv[0],
      nbench);
    printf("Available benchmarks  :\n");

    for(i = 0; i < nbench; i++){
      printf("%i - %s\n", i, benchmarks[i].name);
    }
    for(i = 0; i < nbench2d; i++){
      printf("%i - %s\n", i, benchmarks2d[i].name);
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

    double exec_time = 0.0;
    for(int bs = 0; bs < maskl; bs++){

      if (benchmask[bs] == '1' && bs < nbench) {
        exec_time = 0.0;
        exec_time = test1d(nruns, dimx, dimt, benchmarks[bs]);
      } else if (benchmask[bs] == '1'){
        exec_time = test2d(nruns, dimx, dimy, dimt, benchmarks2d[bs - nbench]);
      }
      printf("TODO : exec_time stats %10.3f", exec_time);
    }
    return 0;
  }

}

double
test2d(int nruns, int dimx, int dimy, int dimt, struct benchspec2d benchmark)
{
  int i;

  double ** data_in = alloc_double_mx(dimx, dimy);
  double ** data_out = alloc_double_mx(dimx, dimy);

  struct args_dimt args = get2dargs(dimx, dimy, dimt, data_in, benchmark);
  struct benchscore scores[nruns];

  for (i = 0; i <= nruns; i ++) {
    if (i == 0) {
      benchmark.variant(args, scores, data_out);
    }
    benchmark.variant(args, scores, data_out);
  }
  // ! TODO
  return 0.0;
}

double
test1d(int nruns, int dimx, int dimt, struct benchspec benchmark)
{
  double * data = (double *) aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(double) * dimx);
  // !! TODO !!!
  return 0.0;
}

struct args_dimt
get2dargs(int dimx, int dimy, int dimt, double ** in_data,
  struct benchspec2d bs)
{
  if(dimx == 0 && dimy == 0 && dimt == 0){
    dimx = bs.width;
    dimy = bs.height;
    dimt = bs.iters;
  }
  struct args_dimt res = {dimx, dimy, dimt,in_data} ;
  return res;
}