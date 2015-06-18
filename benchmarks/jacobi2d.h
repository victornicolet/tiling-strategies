#include "../utils.h"

#ifndef JACOBI2D_H
#define JACOBI2D_H

/* Debugging dimensions */
#ifdef DEBUG
  #define DBG_SIZE 1024
  #define DBG_ITER 512
#endif
// App
//#define SEQ
#define CHECK_ON_SIZE 8
/* TILE DIMENSIONS */
/* Iterations within a tile */
#define T_ITERS 32
/* Fill a cache line with 8 doubles or 16 float */
#define T_WIDTH 32
#define T_HEIGHT 32
#define T_WIDTH_FLT 8
/* Different values for overlapped version */

/* Stencil macros */
#define JACOBI2D_T(t,i,j) ((t)[i][j] + (t)[i][j-1] + (t)[i][j+1] +            \
                          (t)[i+1][j] + (t)[i+1][j-1] + (t)[i+1][j-1] +       \
                          (t)[i-1][j-1] + (t)[i-1][j] + (t)[i-1][j+1]) / 9.0

double djbi2d_seq(struct args_dimt, double **, double **);
double djbi2d_half_diamonds(struct args_dimt, double **, double **);

int check2d_default(int, int, int);

#endif /* JACOBI1D */
