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

// TILE DIMENSIONS -------------------------------
// Iterations within a tile
#define T_ITERS 32
// Fill a cache line with 8 doubles or 16 float
#define T_WIDTH_DBL 8
#define T_WIDTH_FLT 16
// Different values for overlapped version
static int T_WIDTH_DBL_OVERLAP =
 (12*1024 - T_ITERS * T_ITERS)/(sizeof(double)*T_ITERS);
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

void djbi1d_skewed_tiles(int n, int jbi_iters, double ** jbi){
  int r, l, bot, top;

  int stg = (jbi_iters / T_ITERS) + 1;
  int strp = (n / T_WIDTH_DBL); 
  for (long Ti = -1; Ti < strp + 1; Ti++) {
    for (long Tt = 0; Tt < stg; Tt++) {
      // Tile height
      bot = max(Tt * T_ITERS, 1);
      top = min((Tt + 1) * T_ITERS, jbi_iters);
      for(int t = bot; t < top; t++){
        // Line boundaries
        l = max(Ti * T_WIDTH_DBL + (t - bot) , 1);
        r = min((Ti + 1) * T_WIDTH_DBL + (t - bot), n - 1);
        for(int i = l; i < r; i++){
          JBI1D_STENCIL_T(jbi);
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
      int w0 = r0 -l0;
      int l = max((l0 - h ), 1);
      int r = min((r0 + h), n-1);
      int w = r - l;
      {
        // Read tile base
        for(i = l ; i < r ; i++ ){
          tile[0][i-l] = jbi[0][i];
          tile[1][i-l] = 0.0;
        }

        lvl0 = tile[0];
        lvl1 = tile[1];

        for(int t = 0 ; t < h; t++){
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

void djbi1d_skewed_tiles_test(int n, int iters, double ** jbi, 
  struct benchscore * bsc){
  
  double ** jbi_full_mat = allocmatrix_d(iters, n);

  for(int i = 0; i < n; i++){
    jbi_full_mat[0][i] = jbi[0][i];
  }
 printf("r\n");
  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  djbi1d_skewed_tiles(n, iters, jbi_full_mat);
  clock_gettime( CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);

  for(int i = 0; i < n; i++){
    jbi[1][i] = jbi_full_mat[iters-1][i];
  }

  freematrix_d(jbi_full_mat, iters);
}