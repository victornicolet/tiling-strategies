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

double jacobi2d_seq(int iters, int h, int w, double ** img){
  uint8_t x,y,t;

  ALLOC_MX(temp, double, w, h);

  clock_gettime(CLOCK_MONOTONIC, &tbegin);

  for(t = 0; t < iters; t ++){
    for(x = 1; x < w - 1; x ++){
      for(y = 1; y < h - 1; y ++){
        temp[x][y] = JACOBI2D_T(img,x,y);
      }
    }
    // Copy back into the image
    for(x = 0; x < w; x++){
      memcpy(img[x], temp[x], h * sizeof(double));
    }
  }

  clock_gettime(CLOCK_MONOTONIC, &tend);

  FREE_MX(temp, w)

  return ELAPSED_TIME(tend, tbegin)/1000.0 ;
}

int main(int argc, char ** argv){
  return 0;
}