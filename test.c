#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>

#include "jacobi1d.h"
#include "utils.h"

struct benchspec benchmarks[] = {
  {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap_test},
  {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive_test },
  {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test},
  {"JACOBI1D_SWAP_SEQ", djbi1d_swap_seq},
};


#define MIN_POW 10
#define MAX_POW 16

benchmark_jacobi1d(char *);

int main(int argc, char ** argv){

}
