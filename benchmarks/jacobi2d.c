#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include <inttypes.h>
#include "jacobi2d.h"


static struct timespec tend;
static struct timespec tbegin;

/*
 * WIP : need to adjust indexes in computations and load/store
 */
void
djbi2d_half_diamonds(struct args_dimt args, double ** data_in,
 struct benchscore * bsc, double ** data_out)
{
  int iters = args.iters ;
  int stepwidth = 2 * iters ;
  int x_steps = (args.width / stepwidth );
  int y_steps = (args.height / stepwidth );

  /* The tile is stored into temporary arrays : share them ! */
  double ** tile0 = alloc_double_mx(stepwidth, stepwidth);
  double ** tile1 = alloc_double_mx(stepwidth, stepwidth);
  /* Store the layer btw the two stages */
  double *** tmp_layer = malloc(2 * sizeof(*tmp_layer));
  tmp_layer[1] = alloc_double_mx(args.width, args.height);
  tmp_layer[0] = alloc_double_mx(args.width, args.height);

  int xx, yy, x, y, t;
/* First loop : base-down pyramids */
#pragma omp parallel for schedule(static) shared(tile0, tile1, tmp_layer)
  for (xx = 0; xx < x_steps; xx ++) {
    for (yy = 0; yy < y_steps; yy ++) {
/* Tile base :
 *     x0, y1 ---- x1,y1
 *         |           |
 *         |           |
 *       x0,y0 ---- x1,y0
 */
      int x0, x1, y0, y1, x_orig, y_orig;

      x_orig = max(xx * stepwidth, 0);
      y_orig = max(yy * stepwidth, 0);

      for (t = 0; t < iters; t++) {
        x0 = max( xx * stepwidth + t, 0);
        x1 = min((xx + 1) * stepwidth - t, args.width);
        y0 = max( yy * stepwidth + t, 0);
        y1 = min((yy + 1) * stepwidth - t, args.height);

        /* Save results in layer */

        tmp_layer[0][x0][y0] = tile0[x0 - x_orig][y0 - y_orig];
        tmp_layer[0][x1][y1 - 1] = tile0[x1 - x_orig][y1 - y_orig] - 1;

        for (y = y0 + 1; y < y1 - 1; y++) {
          tmp_layer[0][x0][y] = tile0[x0 - x_orig][y - y_orig];
          tmp_layer[1][x0][y] = tile0[x0 - x_orig][y - y_orig];
          tmp_layer[0][x1][y] = tile0[x1 - x_orig][y - y_orig];
          tmp_layer[1][x1][y] = tile0[x1 - x_orig][y - y_orig];
        }

        tmp_layer[0][x0][y0] = tile0[x0 - x_orig][y0 - y_orig];
        tmp_layer[0][x1 - 1][y0] = tile0[x0 - 1 - x_orig][y0 - y_orig];

        for (x = x0 + 1; x < x1 - 1; x++) {
          tmp_layer[0][x][y0] = tile0[x - x_orig][y0 - y_orig];
          tmp_layer[1][x][y0] = tile0[x - x_orig][y0 - y_orig];
          tmp_layer[0][x][y1] = tile0[x - x_orig][y1 - y_orig];
          tmp_layer[1][x][y1] = tile0[x - x_orig][y1 - y_orig];
        }

        /* Stencil computations */
        for (x = x0; x < x1; x++) {
          for (y = y0; y < y1; y++) {
            tile1[x][y] = JACOBI2D_T(tile0, x, y);
          }
        }
        /* Copy result array into source array */
        for (x = x0; x < x1; x++) {
          memcpy(tile0[x], tile1[x], (y1 - y0) * sizeof(*tile0[x]));
        }
      }
    }
  }
/* Second loop : tip-down pyramids */
#pragma omp parallel for schedule(static) shared(tile0, tile1, tmp_layer)
  for (xx = 0; xx < x_steps; xx ++) {
    for (yy = 0; yy < y_steps; yy ++) {
      int x0, x1, y0, y1, x_orig, y_orig;

      x_orig = max(xx * stepwidth + stepwidth / 2, 0);
      y_orig = max(yy * stepwidth + stepwidth / 2, 0);

      for (t = 0; t < iters; t++) {
        x0 = max(x_orig - t, 0);
        x1 = min(x_orig + t, args.width);
        y0 = max(y_orig - t, 0);
        y1 = min(y_orig + t, args.height);

        for (x = x0; x < x1; x++) {
          for (y = y0; y < y1; y++) {
            tile1[x][y] = JACOBI2D_T(tile0, x, y);
          }
        }
      }
    }
  }
  free(tile1);
  free(tile0);
}


void
djbi2d_seq(struct args_dimt args, double ** data_in, struct benchscore * bsc,
  double ** data_out)
{
  uint8_t x, y, t;
  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  //#pragma scop
  for (t = 0; t < args.iters; t ++) {
    for (x = 1; x < args.width - 1; x ++) {
      for (y = 1; y < args.height - 1; y ++) {
        data_out[x][y] = JACOBI2D_T(data_in,x,y);
      }
    }
/* Copy back into the image */
    for (x = 0; x < args.width; x++) {
      memcpy(data_in[x], data_out[x], args.height * sizeof(double));
    }
  }
  //#pragma endscop

  clock_gettime(CLOCK_MONOTONIC, &tend);

  if (bsc != NULL) {
    bsc->wallclock = ELAPSED_TIME(tend, tbegin);
  }
}


int
djbi2d_(struct args_dimt args, double ** input, double ** output)
{
  int i;
  double ** ref_output = alloc_double_mx(args.width, args.height);
  djbi2d_seq(args, input, NULL, ref_output);
  for (i = 0; i < args.width; i ++) {
    if (compare(ref_output[i], output[i], args.height) == 0) {
      free_mx((void **) ref_output, args.width);
      return 0;
    }
  }
  free_mx((void **) ref_output, args.width);
  return 1;
}

int
check2d_default(int dimx, int dimy, int dimt)
{
  return 0;
}