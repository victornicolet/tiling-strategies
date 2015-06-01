#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include "utils.h"

#define UPPER_THRESHOLD 0.8
#define LOWER_THRESHOLD 0.6

#define GAUSS_MASK_5(im,x,y) (2.0*((im)[x-2][y-2]+(im)[x+2][y+2]+             \
    (im)[x+2][y-2]+(im)[x-2][y+2]) +                                          \
    4.0*((im)[x-2][y+1]+(im)[x-1][y+2]+(im)[x+1][y+2]+(im)[x+2][y+1]+         \
      (im)[x+2][y-1]+(im)[x+1][y-2]+(im)[x-1][y-2]+(im)[x-2][y-1])+           \
    9.0*((im)[x-1][y-1]+(im)[x-1][y+1]+(im)[x+1][y+1]+(im)[x+1][y-1])+        \
    12.0*((im)[x][y-1]+(im)[x][y+1]+(im)[x-1][y]+(im)[x+1][y])+               \
    15 * ((im)[x][y])) * (1.0/159.0)

#define GRAD_X(im,x,y) (im)[x+1][y]-(im)[x-1][y]
#define GRAD_Y(im,x,y) (im)[x][y+1]-(im)[x][y-1]
#define GRAD2D(im,x,y) GRAD_Y((im),x,y)+GRAD_X((im),x,y)

#define EDGE(im,x,y) atan((GRAD_X((im),x,y))/(GRAD_Y((im),x,y)))

int
pipeline_canny_naive(int width, int height,
  double ** image, int hit)
{

  uint8_t x,y,t;

  /* Padding image for gaussian filter */
  double ** image_padded;
  alloc_double_mx(image_padded, width, height);

  for (x = 0; x < width; x ++) {
    for (y = 0; y < height; y ++) {
      image_padded[x + 4][y + 4] = image[x][y];
    }
  }


  /* Apply Gaussian filter */
  double ** image_G;
  alloc_double_mx(image_G, width, height);

  for (x = 0; x < width; x ++) {
    for (y = 0; y < height; y ++) {
      image_G[x][y] = GAUSS_MASK_5(image_padded, (x+4), (y+4));
    }
  }
  free_mx(image_padded, width + 4);


  /* Apply gradient intensity */
  double ** gradx, ** grady;
  alloc_double_mx(grady, height, width);
  alloc_double_mx(gradx, height, width);

  for (x = 0; x < width; x ++) {
    for (y = 0; y < height; y ++) {
      grady[x][y] = GRAD_Y(image_G, x, y);
      gradx[x][y] = GRAD_X(image_G, x, y);
    }
  }


  /* Edge detection 1 */
  double ** theta;
  alloc_double_mx(theta, width, height);

  for (x = 0; x < width; x ++) {
    for (y = 0; y < height; y ++) {
      theta[x][y] = EDGE(im, x, y)
    }
  }


  /* Edge levelling */
  double ** tmp;
  alloc_double_mx(tmp, width, height);

  for t = 0; t < hit; t++) {
    for (x = 0; x < width; x ++) {
      for (y = 0; y < height; y ++) {

      }
    }
  }
}

void
gaussian_filter(int width, int height,int x,int y, double ** image)
{

}
