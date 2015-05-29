#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include "utils.h"

#define GAUSS_MASK_5(im,x,y) (2.0*(im[x-2][y-2]+im[x+2][y+2]+\
    im[x+2][y-2]+im[x-2][y+2]) +\
    4.0*(im[x-2][y+1]+im[x-1][y+2]+im[x+1][y+2]+im[x+2][y+1]+\
      im[x+2][y-1]+im[x+1][y-2]+im[x-1][y-2]+im[x-2][y-1])+\
    9.0*(im[x-1][y-1]+im[x-1][y+1]+im[x+1][y+1]+im[x+1][y-1])+\
    12.0*(im[x][y-1]+im[x][y+1]+im[x-1][y]+im[x+1][y])+\
    15 * (im[x][y])) * (1.0/159.0)

#define GRAD_X(im,x,y) im[x+1][y]-im[x-1][y]
#define GRAD_Y(im,x,y) im[x][y+1]-im[x][y-1]
#define GRAD2D(im,x,y) GRAD_Y(im,x,y)+GRAD_X(im,x,y) 

#define LOOP2D(x,y,x0,xn,y0,yn)   for(x = x0; x < xn; x++){ \
    for(y = y0; y < yn; y++){
#define END_LOOP2D } }

#define EDGE(im,x,y) atan((GRAD_X(im,x,y))/(GRAD_Y(im,x,y)))

int pipeline_canny_naive(int width, int height,
  double ** image, int hit){

  uint8_t x,y;

  // Padding image for gaussian filter
  ALLOC_MX(image_padded, width + 4, height + 4)

  LOOP2D(x, y, 0, width, 0, height)
    image_padded[x + 4] = image[x];
  END_LOOP2D
  // Apply Gaussian filter
  ALLOC_MX(image_G, width, height)
  LOOP2D(x, y, 0, width, 0, height)
      image_G[x][y] = GAUSS_MASK_5(image, x+4, y+4);
  END_LOOP2D
  FREE_MX(image_padded, width + 4)
  // Apply gradient intensity
  ALLOC_MX(gradx, width, height)
  ALLOC_MX(grady, width, height)
  LOOP2D(x, y, 0, width, 0, height)
    grady[x][y] = GRAD_Y(image_G, x, y);
    gradx[x][y] = GRAD_X(image_G, x, y);
  END_LOOP2D
  // Edge detection 1
  ALLOC_MX(theta, width, height)
  LOOP2D(x, y, 0, width, 0, height)
    theta[x][y] = EDGE(im, x, y)
  END_LOOP2D
  // Edge levelling
}

void gaussian_filter(int width, int height,int x,int y, double ** image){

}
