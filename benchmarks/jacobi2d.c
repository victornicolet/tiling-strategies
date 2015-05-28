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

int check_correct(double ** input, double ** output){

}

int main(int argc, char ** argv){
  return 0;
}