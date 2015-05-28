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

static int run;

struct jbi2d_args {
  int X;
  int Y;
  int iters;
  double ** input;
};

void djbi2d_half_diamonds(struct jbi2d_args args, struct benchscore * bsc){

  int stepwidth = 2 * args.iters ;
  int x_steps = (args.X / stepwidth );
  int y_steps = (args.Y / stepwidth );

  int xx,yy,x,y,t;
  // First loop : base-down pyramids
  for(xx = 0; xx < x_steps; xx ++){
    for(yy = 0; yy < y_steps; yy ++){
      /* Tile base :
        x0, y1 ---- x1,y1
          |           |
          |           |
        x0,y0 ---- x1,y0
        */
      int x0, x1, y0, y1;

      for(t = 0; t < args.iters; t++){
        x0 = max( xx * stepwidth + t, 0);
        x1 = min((xx + 1) * stepwidth - t, args.X);
        y0 = max( yy * stepwidth + t, 0);
        x1 = min((yy + 1) * stepwidth - t, args.Y);

        for(x = x0; x < x1; x++){
          for(y = y0; y < y1; y++){
            // Do some jacobi
          }
        }
      }
    }
  }

  // Second loop : tip-down pyramids
  for(xx = 0; xx < x_steps; xx ++){
    for(yy = 0; yy < y_steps; yy ++){
      int x0, x1, y0, y1;

      for(t = 0; t < args.iters; t++){
        x0 = max( xx * stepwidth + t, 0);
        x1 = min((xx + 1) * stepwidth - t, args.X);
        y0 = max( yy * stepwidth + t, 0);
        x1 = min((yy + 1) * stepwidth - t, args.Y);

        for(x = x0; x < x1; x++){
          for(y = y0; y < y1; y++){
            // Do some jacobi
          }
        }
      }
    }
  }
}

void djbi2d_seq(struct jbi2d_args args, struct benchscore * bsc){
  uint8_t x,y,t, iters;

  ALLOC_MX(temp, double, args.X, args.Y);

  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  #pragma scop
  for(t = 0; t < iters; t ++){
    for(x = 1; x < args.X - 1; x ++){
      for(y = 1; y < args.Y - 1; y ++){
        temp[x][y] = JACOBI2D_T(args.input,x,y);
      }
    }
    // Copy back into the image
    for(x = 0; x < args.X; x++){
      memcpy(img[x], temp[x], args.Y * sizeof(double));
    }
  }
  #pragma endscop

  clock_gettime(CLOCK_MONOTONIC, &tend);

  FREE_MX(temp, args.X)
}

int djbi2d_(struct jbi2d_args input, double ** output){
  djbi2d_seq(input, NULL);
  for(int i = 0; i < input.X; i++){
    if(compare(input[i], output[i]) == 0){
      return 0;
    }
  }
  return 1;
}

int main(int argc, char ** argv){
  return 0;
}