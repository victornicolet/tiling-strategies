#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>

#include "utils.h"

// Cache line size of 64 bytes on most x86
#define  CACHE_LINE_SIZE 64
#define  L1_CACHE_SIZE 6044
// TILE DIMENSIONS -------------------------------
// Iterations within a tile
#define T_ITERS 32
// Fill a cache line with 8 doubles or 16 float
#define T_WIDTH_DBL 16
#define T_WIDTH_FLT 32
// Different values for overlapped version
static int T_WIDTH_DBL_OVERLAP =
 (L1_CACHE_SIZE - T_ITERS * T_ITERS)/(sizeof(double)*T_ITERS);
 // Specifi for diamond tiles
 #define T_WIDTH_DBL_DIAM 8
 #define T_ITERS_DIAM 32
// -------------------------------
#define JBI1D_STENCIL_T(jbi) jbi[t][i] = \
  (jbi[t-1][i] + jbi[t-1][i-1] + jbi[t-1][i+1]) / 3.0
#define JBI1D_STENCIL_SW(jbi) jbi[1][i] = \
  (jbi[0][i-1] + jbi[0][i] + jbi[0][i]) / 3.0
#define JBI1D_STENCIL(lvl1,lvl0) lvl1[i] = \
  (lvl0[i-1] + lvl0[i+1] + lvl0[i]) / 3.0 

static struct timespec tend;
static struct timespec tbegin;

static int run;


void djbi1d_omp_naive_test(int, int, double**, struct benchscore * );
void djbi1d_omp_overlap_test(int, int, double**, struct benchscore * );
void djbi1d_skewed_tiles_test(int, int, double **, struct benchscore * );

struct benchspec {
  // Name of the benchmark
  char *name;
  // Function to call
  void (*variant)(int, int, double**, struct benchscore *);
};


void djbi1d_skewed_tiles(int strips, int tsteps, double * dashs, \
  double * slashs, int ** task){
    /* In this version we assume T_ITERS >= T_WIDTH_DBL so we have to make a
    difference between regular tiles ( parallelogram-shaped ones) and triangular
    tiles.
    Parallelograms can be cut on one corner (for example leftmost parallelograms
    are cut on top left corner if T_ITERS > T_WIDTH_DBL ).
    With this assumptions parallelograms have always the same dependencies, and
    triangles also.
    */
    #pragma omp parallel
    {
      #pragma omp master
      for(int Ti = 0; Ti < strips; Ti++){
        for(int Tt = 0; Tt < tsteps; Tt++){
          // Strip number
          int strpno = (Ti + Tt * 2) % strips;
          // Bottom tiles : only left-to-right dependencies
          if(Tt == 0){
            #pragma omp task depend(in : task[Ti-1][0]) \
              depend(out : task[Ti + 1][0])
            {
              double l1[T_WIDTH_DBL + 2];
              double l2[T_WIDTH_DBL + 2];
              double * tmp;

              for(int i = 2; i < T_WIDTH_DBL; i++){
                l1[i] = dashs[strpno + i -2];
              }
              for(int t = 0; t < T_ITERS * 2; t+=2){
                l1[0] = slashs[Tt + t];
                l1[1] = slashs[Tt + t + 1];

                for(int i = 2; i < T_WIDTH_DBL + 2; i++){
                  l2[i] = (l1[i -2] + l1[i - 1] + l1[i]) / 3.0 ;
                }

                slashs[Tt + t] = l2[T_WIDTH_DBL];
                slashs[Tt + t + 1] = l2[T_WIDTH_DBL + 1];

                *tmp = *l1;
                *l1 = *l2;
                *l2 = *tmp;
              }
            }

          } else if(Ti == 0){
            /* Left edge tile : triangular
            Only one in dependency, one out
            Here assume T_ITERS > T_WIDTH_DBL
            */
            #pragma omp task depend(in : task[Ti][Tt-1]) \
              depend(out : task[Ti+1][Tt])
            {
              double l1[T_WIDTH_DBL + 2];
              double l2[T_WIDTH_DBL + 2];
              double *tmp; 

              for(int i = 1; i < T_WIDTH_DBL + 2; i++){
               l1[i] = dashs[strpno + i - 1];
              }

              for(int t= 0; t < 2 * T_ITERS; t+=2){
                int right = max(T_WIDTH_DBL - t, 0);
                l1[0] = 0;
                for(int i = 1; i < right; i++){
                  l2[i] = (l1[i - 1] + l1[i] + l1[i + 1]) / 3.0;
                }
                slashs[Tt + t] = l2[right];
                slashs[Tt + t + 1] = l2[right - 1];

                *tmp = *l1;
                *l1 = *l2;
                *l2 = *tmp; 
              }
              for(int i = 2; i < T_WIDTH_DBL; i++){
                dashs[strpno + i - 2] = l1[i];
              }

              
            }
          } else {
            // Regular tile
            // Two in and out dependencies
            #pragma omp task depend(in: task[Ti-1][Tt], task[Ti][Tt-1]) \
              depend(out : task[Ti + 1][Tt], task[Ti][Tt + 1])
            {
              double l1[T_WIDTH_DBL+2];
              double l2[T_WIDTH_DBL+2];
              double * tmp;

              // Load dash in the stack
              for(int i = 2; i < T_WIDTH_DBL + 2; i++){
                l1[i] = dashs[strpno + i - 2];
              }
              for(int t = 0; t < T_ITERS * 2; t+=2){
                // Load slash part
                l1[0] = slashs[Tt + t];
                l1[1] = slashs[Tt + t + 1];
                for(int i = 2; i < T_WIDTH_DBL + 2; i++){
                  l2[i] = (l1[i -2] + l1[i - 1] + l1[i]) / 3.0 ;
                }
                // Write slash
                slashs[Tt + t] = l2[T_WIDTH_DBL];
                slashs[Tt + t + 1] = l2[T_WIDTH_DBL + 1];

                *tmp = *l1;
                *l1 = *l2;
                *l2 = *tmp;
              }
              // Write dash
              for(int i = 2; i < T_WIDTH_DBL; i++){
                dashs[strpno + i - 2] = l1[i];
              }
            }
          }
        }
      }
  }
}



void djbi1d_diamond_tiles(int n,int jbi_iters, double ** jbi, 
  struct benchscore * bsc){
  int r, l, bot, top;

  int stg = (jbi_iters / T_ITERS_DIAM) + 1;
  int strp = (n / T_WIDTH_DBL_DIAM); 
  for (long Ti = -1; Ti < strp + 1; Ti++) {
    for (long Tt = 0; Tt < stg; Tt++) {
      // Tile height
      bot = max(Tt * T_ITERS_DIAM, 1);
      top = min((Tt + 1) * T_ITERS_DIAM, jbi_iters);
      for(int t = bot; t < top; t++){
        // Line boundaries
        l = max(Ti * T_WIDTH_DBL_DIAM - (t - bot) , 1);
        r = min((Ti + 1) * T_WIDTH_DBL_DIAM - (t - bot), n - 1);
        for(int i = l; i < r; i++){
          JBI1D_STENCIL_T(jbi);
        }
      }
    }
  }
}

void djbi1d_omp_naive(int n, int jbi_iters, double ** jbi){

  int t,i;
  double * tmp;
  for(t = 0; t < jbi_iters-1; t++){
    #pragma omp parallel for schedule(static)
    for(i = 1; i < n - 1; i++){
      #pragma  ivdep
      JBI1D_STENCIL_SW(jbi);
    }
    tmp = jbi[0];
    jbi[0] = jbi[1];
    jbi[1] = tmp;
  }
}

void djbi1d_omp_overlap(int n, int jbi_iters, int stages, double ** jbi){
  int Ti,Tt,t,i;
  for( Tt = 0; Tt < stages - 1; Tt ++){
    double * tmp;
    #pragma omp parallel for schedule(static)
    for(Ti = 0; Ti <= (n / T_WIDTH_DBL_OVERLAP); Ti ++){
      double tile[2][T_WIDTH_DBL_OVERLAP + 2 * T_ITERS];
      double * lvl0, *lvl1, *tmp;

      // Compute tile bounds
      int bot = max((T_ITERS * Tt), 0);
      int top = min((T_ITERS * (Tt + 1)), jbi_iters);
      int h = top - bot;
      int l0 = max((T_WIDTH_DBL_OVERLAP * Ti), 1);
      int r0 = min((T_WIDTH_DBL_OVERLAP * (Ti + 1)), n-1);
      int l = max((l0 - h ), 1);
      int r = min((r0 + h), n-1);
      {
        // Read tile base
        for(i = l ; i < r ; i++ ){
          tile[0][i-l] = jbi[0][i];
          tile[1][i-l] = 0.0;
        }

        lvl0 = tile[0];
        lvl1 = tile[1];

        for(t = 0 ; t < h; t++){
          int lt = max(t , 0);
          int rt = min((r - l - t), n);
          for(int i = lt ; i < rt; i++){
            JBI1D_STENCIL(lvl1,lvl0);
          }
          tmp = lvl0;
          lvl0 = lvl1;
          lvl1 = tmp;
        }

        // Write tile top
        for(i = l0 ; i < r0 ; i++ ){
          jbi[1][i] = lvl1[i-l];
        }
      }
    }
    tmp = jbi[0];
    jbi[0] = jbi[1];
    jbi[1] = tmp;
  }

}


void djbi1d_swap_seq(int n, int jbi_iters, double ** jbi, 
  struct benchscore * bsc){
  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  int t,i;
  double *a, *b, *tmp;
  a = jbi[0];
  b = jbi[1];
  for(t = 1; t < jbi_iters; t++){
    for(i = 1; i < n - 1; i++){
      JBI1D_STENCIL(b, a);
    }
    tmp = b;
    b = a;
    a = tmp;
  }
  jbi[1] = a;
  // End
  clock_gettime( CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
}



struct benchspec benchmarks[] = {
  {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap_test},
  {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive_test },
  {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test},
  {"JACOBI1D_SWAP_SEQ", djbi1d_swap_seq},
};



int main(int argc, char ** argv){

    int nbench = sizeof(benchmarks) / sizeof(struct benchspec);

    if(strcmp(argv[1], "help") == 0){
        printf("Usage: %s <Nruns> <Mask : %i> [ <Width> <Time iterations>]\n", 
        argv[0], nbench);
        printf("Build mask with : OVERLAP NAIVE SKEWED_TILES SEQUENTIAL\n");
        printf("Dimensions : \n");
        printf("T_ITERS : \t\t%i\nT_WIDTH_DBL : \t\t %i\n", T_ITERS, 
          T_WIDTH_DBL);
        printf("T_WIDTH_DBL_OVERLAP : \t%i\n", T_WIDTH_DBL_OVERLAP);
        printf("T_WIDTH_DBL_DIAM : \t%i\n", T_WIDTH_DBL_DIAM );
        return 0;
    }

    if(argc < 3){
      printf("Usage: %s <Nruns> <Mask : %i> [ <Width> <Time iterations>]\n", 
        argv[0], nbench);
      return 0;
    }

    int i, j, iter;

    int tab_size, jbi_size;
    if(argc == 5){
      tab_size = atoi(argv[3]);
      jbi_size = atoi(argv[4]);
    } else {
      tab_size = 16*4096;
      jbi_size = 1024;
    }
    tab_size = ((tab_size - 1) / (2*T_ITERS - 1)) * (2*T_ITERS - 1);


    double ** jbi = (double **) malloc(sizeof(double) * 2);
    for(i = 0; i < 2; i++){
      jbi[i] = (double *) malloc(sizeof(double) * tab_size);
    }
    for(j = 0; j < tab_size; j++){
      jbi[0][j] = (tab_size  - j)* j / 100.0  ;
    }

    char *benchmask = argv[2];
    if(strlen(benchmask) != nbench){
      printf("Error : not a valid mask ! Your mask must be %i bits long\n", 
        nbench);
      return -1;
    }

    int nruns = atoi(argv[1]);

    printf("Input : \n");
    for(i = 0; i < min(tab_size,8); i++){
      printf("%10.3f", jbi[0][i]);
    }

    printf("\n");
    double accu;
    for(int bs = 0; bs < nbench; bs++){
      if (benchmask[bs] == '1') {
        struct benchscore score[nruns + 1];
        accu = 0.0;
        for(iter = 0; iter < nruns + 1; iter++){
          score[iter].name = benchmarks[bs].name;
          benchmarks[bs].variant(tab_size, jbi_size, jbi, &score[iter]);
          
          if(iter > 0) {
            printf("%s : Run %i ...", score[iter].name, iter + 1 );
            printf("\t\t %13f ms\n", score[iter].wallclock * 1000.0 );
            accu += score[iter].wallclock ;
          }
        }
        printf("\n------------- %s ---------\n", benchmarks[bs].name);
        printf("Result: \n");
        for(i = 0; i < min(tab_size ,8); i++){
          printf("%10.3f", jbi[1][i]);
        }
        printf("\n----------------------\n");
        printf("Total time :\t %13f ms\n", (double) accu * 1000.0);
        printf("Average time :\t %13f ms\n", 
          (double) (accu * 1000.0 / (nruns)));

        // Reinitialisation de la matrice
        for(j = 0; j < tab_size; j++){
          jbi[0][j] = (tab_size  - j)* j / 100.0  ;
        }
      }
    }

    free(jbi[1]);
    free(jbi[0]);
    free(jbi);

}


void djbi1d_omp_naive_test(int n, int jbi_iters, double ** jbi,
  struct benchscore * bsc){

  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  djbi1d_omp_naive(n ,jbi_iters, jbi);
  clock_gettime(CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
  for(int i = 0; i < n; i++){
    jbi[1][i] = jbi[0][i];
  }
}

// For overlap, return the time spent running jbi - memory allocation no included
void djbi1d_omp_overlap_test(int n, int iters, double ** jbi, 
  struct benchscore * bsc){
  int stages = (iters / T_ITERS) + 1;

  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  djbi1d_omp_overlap(n, iters, stages, jbi);
  clock_gettime( CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);

  for(int i = 0; i < n; i++){
    jbi[1][i] = jbi[0][i];
  }
}

void djbi1d_skewed_tiles_test(int n, int iters, double ** jbi, \
  struct benchscore * bsc){
  
  int timesteps = (iters / T_ITERS) + 1;
  int strips = (n / T_WIDTH_DBL);
  double * jbi_dashs = (double*) malloc(sizeof(double) * strips * T_WIDTH_DBL);
  double * jbi_slashs = (double*) malloc(sizeof(double) * timesteps * T_ITERS);

  ALLOC_MX(tasks, int, strips, timesteps)

  for(int i =0; i < n; i++){
    jbi_dashs[i] = jbi[0][i];
  }

  clock_gettime( CLOCK_MONOTONIC, &tbegin);

  djbi1d_skewed_tiles(strips , timesteps, jbi_dashs, jbi_slashs, tasks);

  clock_gettime( CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);

  FREE_MX(tasks, strips)

  free(jbi_dashs);
  free(jbi_slashs);
}