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

struct benchspec benchmarks[] = {
  {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap},
  {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive },
  {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test},
  {"JACOBI1D_SK_FULL_TILES", djbi1d_sk_full_tiles_test},
  {"JACOBI1D_SWAP_SEQ", djbi1d_swap_seq},
  {"JACOBI2D_SEQ", djbi2d_seq}
};


#define MIN_POW 10
#define MAX_POW 16

int main(int argc, char ** argv){

  int nbench = sizeof(benchmarks) / sizeof(struct benchspec);

  if(argc < 3){
    printf("Usage : %s <Test runs / alg.> <Mask (%i bit wide) > \n", argv[1], 
      nbench);
    return 0;
  }

}
