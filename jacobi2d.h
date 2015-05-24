#include "utils.h"

#ifndef JACOBI2D
#define JACOBI2D

// Debugging
#define DEBUG

#ifdef DEBUG
  #define DBG_SIZE 1024
  #define DBG_ITER 512
 #endif
// App
#define SEQ
#define CHECK_ON_SIZE 8
// TILE DIMENSIONS -------------------------------
// Iterations within a tile
#define T_ITERS 32
// Fill a cache line with 8 doubles or 16 float
#define T_WIDTH_DBL 32
#define T_HEIGHT_DBL 32
#define T_WIDTH_FLT 8
// Different values for overlapped version

#define JACOBI2D_T(t,i,j) (t[i][j] + t[i][j-1] + t[i][j+1] + t[i+1][j] + \
  t[i+1][j-1] + t[i+1][j-1] + t[i-1][j-1] + t[i-1][j] + t[i-1][j+1])/9.0


void djbi1d_omp_naive_test(int, int, double**, struct benchscore * );
void djbi1d_omp_overlap_test(int, int, double**, struct benchscore * );
void djbi1d_skewed_tiles_test(int, int, double **, struct benchscore * );



#endif /* JACOBI1D */