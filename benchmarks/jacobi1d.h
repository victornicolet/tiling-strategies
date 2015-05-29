#include "../utils.h"

#ifndef JACOBI1D_H
#define JACOBI1D_H
// Debugging
//#define DEBUG
//#define DEBUG_GDB
// Check dependencies
//define DEBUG_PARALLEL
// Small dimensions for debugging
#ifdef DEBUG
  #define DBG_SIZE 4096
  #define DBG_ITER 4
#endif
// App
// SEQ -> compile sequential version of the algorithm
// All 'pragmas omp' are removed
//#define SEQ
// Number of doubles displayed
#define DISPLAY_SIZE 8
// TILE DIMENSIONS -------------------------------
// Iterations within a tile
#define T_ITERS 32
// Default tile size
#define T_WIDTH_DBL 32
//#define T_WIDTH_FLT 8
// Different values for overlapped version
static int T_WIDTH_DBL_OVERLAP = T_WIDTH_DBL;
 //(L1_CACHE_SIZE - T_ITERS * T_ITERS)/(sizeof(double)*T_ITERS);
 // Specific for diamond tiles
 #define T_WIDTH_DBL_DIAM 8
 #define T_ITERS_DIAM 32
// -------------------------------
#define JBI1D_STENCIL_T(jbi) jbi[t][i] = \
  (jbi[t-1][i] + jbi[t-1][i-1] + jbi[t-1][i+1]) / 3.0
#define JBI1D_STENCIL_SW(jbi) jbi[1][i] = \
  (jbi[0][i-1] + jbi[0][i] + jbi[0][i]) / 3.0
#define JBI1D_STENCIL(lvl1,lvl0) lvl1[i] = \
  (lvl0[i-1] + lvl0[i+1] + lvl0[i]) / 3.0
//-----------------------------------


void djbi1d_omp_naive(int, int, double**, struct benchscore * );
void djbi1d_omp_overlap(int, int, double**, struct benchscore * );
void djbi1d_skewed_tiles_test(int, int, double **, struct benchscore * );
void djbi1d_sk_full_tiles_test(int, int, double **,struct benchscore * );
void djbi1d_half_diamonds_test(int, int, double **,struct benchscore * );
void djbi1d_sequential(int, int, double **,struct benchscore *);

int check_low_iter(int, int);
int check_tilable(int, int);
int check_default(int, int);

double *
alloc_line(int num_elements)
{
  double * line = aligned_alloc(CACHE_LINE_SIZE, num_elements * sizeof(*line));
  if (line == NULL) {
    fprintf(stderr, "Error while allocating line of %i doubles.\n",
            num_elements);
  }
  return line;
}

#endif /* JACOBI1D_H */
