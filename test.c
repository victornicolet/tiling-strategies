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

void test1d(struct benchspec, int, int);
void test2d(struct benchspec, int, int, int);

int main(int argc, char ** argv){

  int i,t;

  struct benchspec benchmarks[] = {
    {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap, 1},
    {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive, 1},
    {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test, 1},
    {"JACOBI1D_SK_FULL_TILES", djbi1d_sk_full_tiles_test, 1},
    {"JACOBI1D_SWAP_SEQ", djbi1d_swap_seq, 1},
    {"JACOBI2D_SEQ", djbi2d_seq, 2}
  };

  int nbench = sizeof(benchmarks) / sizeof(struct benchspec);

  /*---------- Parameters section -----------*/

  if(argc < 2){
    printf("Usage : %s <Mask (%i bit wide) > <Test runs / alg.>\n", argv[0],
      nbench);
    printf("Available benchmarks  :\n");

    for(i = 0; i < nbench; i++){
      printf("%i - %s\n", i, benchmarks[i].name);
    }

    char * benchmask = argv[1];
    int nruns;
    int maskl;

    if(argc == 3){
      nruns = atoi(argv[2]);
    }
    if( maskl = strlen(benchmask) > nbench){
      printf("Error : not a valid mask ! Your mask must be %i bits long\n",
        nbench);
      return -1;
    }

    int dimx, dimy, dimt;
    if(argc == 4){
      dimx = atoi(argv[3]);
    }
    if(argc == 5){
      dimy = atoi(argv[4])
    }

    /* -------------------------------------- */

    for(int bs = 0; bs < maskl; bs++){
      if (benchmask[bs] == '1') {
        int dim = benchmarks[bs].dim;
        if(benchmarks[bs].dim == 1){
          test1d(benchmarks[bs], dim1, dimt);
        } else if(benchmarks[bs].dim == 2){
          test2d(benchmarks[bs], dim1, dim2, dimt);
        }
      }
    }

    return 0;
  }

}

void test2d(struct benchspec benchmarks, int dimx, int dimy, int dimt){
  ALLOC_MX(data, double, dimx, dimy)
}

void test1d(struct benchspec benchmarks, int dimx, int dimt){
  double * data = (double *) aligned_alloc(CACHE_LINE_SIZE, \
    sizeof(double) * dimx);

}
