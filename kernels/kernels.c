#include <stdlib.h>

// Dimensions
#define SMALL (1 << 5)
#define DEFAULT ( 1 << 10)
#define BIG (1 << 13)
#define HUGE (1 << 16)

int main(int argc, char ** argv){
  return 0;
}

void jacobi1d_kernel(){
  int i,t;
  static int t_max, i_max;

  t_max = DEFAULT;
  i_max = BIG;

  double * m1[i_max];
  double * m0[i_max];

#pragma scop
  for(t = 0; t < t_max; t++){
    for(i = 0; i < i_max; i++){
      m1[i] = (m0[i - 1] + m0[i] + m0[i + 1]) / 3.0;
    }
    //swap elements of the arrays
    for(i = 0; i < i_max; i ++){
      m0[i] = m1[i];
    }
  }
#pragma endscop
}

void jacobi2d_kernel(){
  int i, j, k, l, t;
  int t_max, i_max, j_max;

  t_max = DEFAULT;
  i_max = BIG;
  j_max = BIG;

  double ** m1;
  double ** m0;

#pragma scop
  for(t = 0; t < t_max; t++){
    for(i = 1; i < i_max - 1; i++){
      for(i = 1; j < j_max - 1; j++){
        m1[i][j]= 0.2*(m0[i][j]+m0[i][j-1]+m0[i][1+j]+m0[1+i][j]+m0[i-1][j]);
      }
    }
    //swap elements of the arrays
    for(l = 1; l < i_max - 1; l ++){
      for(k= 1; k < j_max - 1; k++){
        m0[l][k] = m1[l][k];
      }
    }
  }
#pragma endscop
}