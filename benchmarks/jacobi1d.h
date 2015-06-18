#include "../utils.h"

#ifndef JACOBI1D_H
#define JACOBI1D_H

/* Macros used for testing :
 * DEBUG
 * DEBUG_PARALLEL --> check dependencies
 * DEBUG_GDB
 */


/* SEQ -> compile sequential version of the algorithm
 * All 'pragmas omp' are removed
 */

/* TILE DIMENSIONS
 * Iterations within a tile
 */
#define T_ITERS 32
/* Default tile size */
#define T_WIDTH_DBL 32
//#define T_WIDTH_FLT 8
/* Different value for overlapped version */
#define T_WIDTH_DBL_OVERLAP (T_WIDTH_DBL)
 //(L1_CACHE_SIZE - T_ITERS * T_ITERS)/(sizeof(double)*T_ITERS);

/* Specific for diamond tiles */
#define T_WIDTH_DBL_DIAM 8
#define T_ITERS_DIAM 32

#ifndef GROUP_FACTOR
#define GROUP_FACTOR 1
#endif

/* ------------------------------- */
/* Stencil macros */

#define JBI1D_STENCIL_T(jbi) jbi[t][i] = \
  (jbi[t-1][i] + jbi[t-1][i-1] + jbi[t-1][i+1]) / 3.0
#define JBI1D_STENCIL_SW(jbi) jbi[1][i] = \
  (jbi[0][i-1] + jbi[0][i] + jbi[0][i]) / 3.0
#define JBI1D_STENCIL(lvl1,lvl0) lvl1[i] = \
  (lvl0[i-1] + lvl0[i+1] + lvl0[i]) / 3.0
/*-----------------------------------*/


double djbi1d_omp_naive(struct args_dimt, double *, double *);
double djbi1d_omp_overlap(struct args_dimt, double *, double *);
double djbi1d_skewed_tiles_test(struct args_dimt, double *,double *);
double djbi1d_sk_full_tiles_test(struct args_dimt, double *,double *);
double djbi1d_half_diamonds_test(struct args_dimt, double *,double *);
double ljbi1d_half_diamonds_test(struct args_dimt, long *,long *);
double djbi1d_sequential(struct args_dimt, double *,double *);
double ljbi1d_sequential(struct args_dimt, long *,long *);
double djbi1d_hdiam_grouped_test(struct args_dimt, double *, double *);
double djbi1d_hdiam_tasked_test(struct args_dimt, double *, double *);

int check_low_iter(int, int);
int check_tilable(int, int);
int check_default(int, int);

#endif /* JACOBI1D_H */
