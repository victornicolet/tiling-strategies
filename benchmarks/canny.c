#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>

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

#define EDGE(im,x,y) atan((GRAD_X(im,x,y))/(GRAD_Y(im,x,y)))

int pipelline_canny(int width, int height,
  double ** image, int hit){
  // Apply Gaussian filter
  // Apply gradient intesity
  // Edge detection
  // Non-maxima suppression
  // Edge levelling
}

void gaussian_filter(int width, int height,int x,int y, double ** image){

}
