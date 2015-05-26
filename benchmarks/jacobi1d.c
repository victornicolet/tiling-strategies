#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include <inttypes.h>
#include "jacobi1d.h"

static struct timespec tend;
static struct timespec tbegin;

static FILE * csv_file;

static int run;

/* Functions describing different tasks in the computation
  Always inlined in the main body */
void do_i_t(double * dashs, double* slashs, int strpno, int Tt) 
  __attribute__((always_inline));
void do_i0_t(double * dashs, double* slashs, int strpno, int Tt)
  __attribute__((always_inline));
void do_i0_t0(double * dashs, double* slashs, int T, int I)
  __attribute__((always_inline));
void do_in_t(double * dashs, double* slashs, int strpno, int Tt)
  __attribute__((always_inline));


uint8_t ** djbi1d_sk_full_tiles(int strips, int steps, double * dashs, \
  double * slashs){
  /* In this algorithm we work only on full parallelograms : no partial tiles. 
  We are interested in performance measures, and understanding, rather 
  than corectness here */

  ALLOC_MX(tasks, uint8_t, steps, strips)

#ifndef SEQ
  #pragma omp parallel
#endif
  {
#ifndef SEQ
    #pragma omp single
#endif
    {
      int i, t;

      for(t = 0; t < steps; t ++){
        for(i = 1; i < strips; i++){
          tasks[t][i] = 0;
          int strpno = i + t;
          // Stratup task
          if(t == 0 && i == 1){
#ifndef SEQ
            #pragma omp task firstprivate(i,t) \
            depend(out : tasks[t][i])
#endif
            {
              do_i_t(dashs, slashs, 1, 0);
              tasks[t][i] ^= 1;
            }
          }else if(t == 0 && i > 1){
#ifndef SEQ
            #pragma omp task firstprivate(i,t) \
            depend(out : tasks[t][i])\
            depend(in : tasks[i-1][0])
#endif
            {
              do_i_t(dashs, slashs, strpno, 0);
              tasks[t][i] ^= 1;
            }
          } else {
#ifndef SEQ            
            #pragma omp task firstprivate(i,t) \
            depend(out : tasks[t][i]) \
            depend(in  : tasks[t][i-1], tasks[t-1][i+1])
 #endif            
            {
              do_i_t(dashs, slashs, strpno, t);
              tasks[t][i] ^= 1;
            }
          }
        }
      }
    }
  }

  return tasks;
}

void djbi1d_skewed_tiles(int strips, int steps, double * dashs, \
  double * slashs){
    /* In this version we assume T_ITERS = T_WIDTH_DBL so we have to make a
    difference between regular tiles ( parallelogram-shaped ones) and triangular
    tiles, but the pattern is quite regular and dependencies are straightforward
    at the boundaries of the domain.
    ---------------

    steps = (iters / T_ITERS) + 1;
    strips = (n / T_WIDTH_DBL) + 1 ;
    */
    uint8_t tasks[strips+1][steps];
    int Ti, Tt;
#ifndef SEQ
  #pragma omp parallel
#endif
    {
#ifndef SEQ
  #pragma omp master
#endif

      for(Tt = 0; Tt < steps; Tt++){
        for(Ti = 0; Ti < strips; Ti++){
          // Strip number
          int strpno = (Ti + Tt);
          //int strpno = strpno % strips;
          // Slash beginning 
          //int s_index = Tt * T_ITERS * 2;
          // Dash beginning
          //int d_index = strpno * T_WIDTH_DBL;
          // Initial tile
          if( Tt == 0 && Ti == 0){
#ifndef SEQ
  #pragma omp task firstprivate(Ti,Tt) \
    depend(out : tasks[0][0])
#endif
            {
              do_i0_t0(dashs, slashs, Ti, Tt);
            }
          } else if(Tt == 0 && Ti < strips - 1){
            // Bottom tiles : only left-to-right dependencies + top out
#ifndef SEQ
  #pragma omp task firstprivate(Ti,Tt) \
    depend(in : tasks[Ti-1][0]) \
    depend(out : tasks[Ti][Tt])
#endif
            {
              do_i_t(dashs, slashs, strpno, 0);
            }

          } else if(Ti == 0 && Tt > 0){
            /* Left edge tile : triangular
            Only one in dependency, one out
            Here assume T_ITERS > T_WIDTH_DBL
            */
#ifndef SEQ
  #pragma omp task firstprivate(Ti,Tt) \
    depend(in : tasks[1][Tt-1]) \
    depend(out : tasks[Ti][Tt])
#endif
            { 
              do_i0_t(dashs, slashs, strpno, Tt);
            }
          } else if(Ti == strips - 1){
#ifndef SEQ
  #pragma omp task firstprivate(Ti,Tt) \
    depend(in: tasks[strips][Tt]) \
    depend(out : tasks[Ti][Tt])
#endif
            {
              do_in_t(dashs, slashs, strpno, Tt);
            }

          }else{
            // Regular tile
            // Two in and out dependencies
#ifndef SEQ
  #pragma omp task firstprivate(Ti,Tt) \
    depend(in : tasks[Ti-1][Tt],tasks[Ti+1][Tt-1]) \
    depend(out : tasks[Ti][Tt])
#endif
            {
              do_i_t(dashs, slashs, strpno, Tt);
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

void djbi1d_omp_naive(int n, int jbi_iters, double ** jbi,
  struct benchscore * bsc){

  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  int t,i;
  double * l1 = (double *) aligned_alloc(CACHE_LINE_SIZE, sizeof(double) * n);
  double * l2 = (double *) aligned_alloc(CACHE_LINE_SIZE, sizeof(double) * n);
  memcpy(l1, jbi[0], n * sizeof(double));

  clock_gettime( CLOCK_MONOTONIC, &tbegin);

  for(t = 1; t < jbi_iters; t++){
#ifdef SEQ
   #pragma omp parallel for schedule(static)
#endif
    for(i = 1; i < n - 1; i++){
      JBI1D_STENCIL(l2, l1);
    }
    memcpy(l1, l2, n * sizeof(double));
  }
  // End
  clock_gettime( CLOCK_MONOTONIC, &tend);

  for(int i = 0; i < n; i++){
    jbi[1][i] = l1[i];
  }

  free(l1);
  free(l2);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
}

void djbi1d_omp_overlap(int n, int jbi_iters, double ** jbi,
  struct benchscore * bsc){
  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  int Ti,Tt,t,i;
  double * tmp;

  for( Tt = 0; Tt <= jbi_iters/T_ITERS; Tt ++){
#ifndef SEQ
    #pragma omp parallel for schedule(static)
#endif
    for(Ti = 0; Ti <= (n / T_WIDTH_DBL_OVERLAP); Ti ++){

      double * lvl1 = (double *) malloc(sizeof(double) * (T_WIDTH_DBL_OVERLAP+
          T_ITERS * 2));
      double * lvl0 = (double *) malloc(sizeof(double) * (T_WIDTH_DBL_OVERLAP+
          T_ITERS * 2));

      double *tmp;

      // Compute tile bounds
      int bot = max((T_ITERS * Tt), 0);
      int top = min((T_ITERS * (Tt + 1)), jbi_iters);
      int h = top - bot;
      int l0 = max((T_WIDTH_DBL_OVERLAP * Ti), 0);
      int r0 = min((T_WIDTH_DBL_OVERLAP * (Ti + 1)), n);
      int l = max((l0 - h ), 0);
      int r = min((r0 + h), n);
      int w = r - l;
      {
        // Read tile base
        for(i = l ; i < r ; i++ ){
          lvl0[i- l ] = jbi[0][i];
          lvl1[i- l] = 0.0;
        }

        for(t = 0 ; t < h; t++){
          int lt = max(t + 1 , 1);
          int rt = min((w - t - 1), (T_WIDTH_DBL_OVERLAP+
          T_ITERS * 2)); 
          for(i = lt ; i < rt; i++){
            JBI1D_STENCIL(lvl1,lvl0);
          }
          memcpy(lvl0, lvl1, 
            (T_WIDTH_DBL_OVERLAP + T_ITERS * 2) * sizeof(double) );

          #ifdef DEBUG
          if(Ti == 0){
              fprintf(csv_file, "Left column ; %i ;%i", Tt, t);
              for( i = 0; i < 8 ; i++){
                fprintf(csv_file, ";%10.3f", lvl1[i] - jbi[0][i + l]);
              }
              fprintf(csv_file, "\n");
          }
          #endif
        }

        // Write tile top
        for(i = l0 ; i < r0 ; i++ ){
          jbi[1][i] = lvl0[i-l];
        }
        free(lvl1);
        free(lvl0);
      }
    }
    memcpy(jbi[0], jbi[1], n * sizeof(double));
  }
  clock_gettime( CLOCK_MONOTONIC, &tend);


  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
}


void djbi1d_swap_seq(int n, int jbi_iters, double ** jbi, 
  struct benchscore * bsc){
    // Boundaries initial condition
  int t,i;
  double * l1 = (double *) aligned_alloc(CACHE_LINE_SIZE, sizeof(double) * n);
  double * l2 = (double *) aligned_alloc(CACHE_LINE_SIZE, sizeof(double) * n);
  memcpy(l1, jbi[0], n * sizeof(double));

  clock_gettime( CLOCK_MONOTONIC, &tbegin);

  for(t = 1; t < jbi_iters; t++){
    for(i = 1; i < n - 1; i++){
      JBI1D_STENCIL(l2, l1);
    }
    memcpy(l1, l2, n * sizeof(double));
  }
  // End
  clock_gettime( CLOCK_MONOTONIC, &tend);
  for(int i = 0; i < n; i++){
    jbi[1][i] = l1[i];
  }

  free(l1);
  free(l2);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
}


struct benchspec benchmarks[] = {
  {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap},
  {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive },
  {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test},
  {"JACOBI1D_SK_FULL_TILES", djbi1d_sk_full_tiles_test},
  {"JACOBI1D_SWAP_SEQ", djbi1d_swap_seq},
};



int main(int argc, char ** argv){

    int nbench = sizeof(benchmarks) / sizeof(struct benchspec);

    if(argc < 3){
      printf("Usage: %s <Nruns> <Mask : %i> [ <Width> <Time iterations>]\n", 
        argv[0], nbench);
      return 0;
    }

    csv_file = fopen("jacobi1d.csv", "w");
    if (csv_file == NULL) {
      exit(1);
    }

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

    int i, j, iter;
    int nruns = atoi(argv[1]);
    int tab_size, jbi_size;
    if(argc == 5){
      tab_size = 1 << atoi(argv[3]);
      jbi_size = 1 << atoi(argv[4]);
    } else {
      #ifdef DEBUG
        tab_size = DBG_SIZE;
        jbi_size = DBG_ITER;
      #else
        tab_size = 1 << 16;
        jbi_size = 1 << 10;
      #endif
    }

    #ifndef DEBUG
      tab_size = ((tab_size - 1) / (2*T_ITERS - 1)) * (2*T_ITERS - 1);
    #endif

    double ** jbi = (double **) malloc(sizeof(double) * 2);
    for(i = 0; i < 2; i++){
      jbi[i] = (double *) malloc(sizeof(double) * tab_size);
    }

    JBI_INIT(jbi, tab_size)

    char *benchmask = argv[2];
    if(strlen(benchmask) != nbench){
      printf("Error : not a valid mask ! Your mask must be %i bits long\n", 
        nbench);
      return -1;
    }

    printf("Input : \n");
    for(i = 0; i < CHECK_ON_SIZE; i++){
      printf("%10.3f", jbi[0][i]);
    }

    #ifdef DEBUG
      fprintf(csv_file, "Whare ?; Data\n");
      fprintf(csv_file, "Input;\n");
      for(i = 0; i < 1 << 8; i++){
        fprintf(csv_file, "%10.3f ;", jbi[0][i]);
      }
      fprintf(csv_file, "\n");
    #endif

    // Get the correct result
    double * check_res = (double *) malloc(sizeof(double) * tab_size);
    struct benchscore bsc;
    djbi1d_swap_seq(tab_size, jbi_size, jbi, &bsc);

    for(int i = 0; i < tab_size; i++){
      check_res[i] = jbi[1][i];
    }

    printf("\n");
    double accu;

    for(int bs = 0; bs < nbench; bs++){
      if (benchmask[bs] == '1') {
        #ifdef SEQ
          printf("WARNING : SEQ defined\n");
        #endif
        struct benchscore score[nruns + 1];
        accu = 0.0;
        for(iter = 0; iter < nruns + 1; iter++){

          JBI_INIT(jbi, tab_size)
          score[iter].name = benchmarks[bs].name;
          benchmarks[bs].variant(tab_size, jbi_size, jbi, &score[iter]);
          
          if(iter > 0) {
            printf("%s : Run %i ...", score[iter].name, iter );
            printf("\t\t %13f ms\n", score[iter].wallclock * 1000.0 );
            accu += score[iter].wallclock ;
          }
        }
        printf("\n------------- %s ---------\n", benchmarks[bs].name);
        printf("Result: \n");
        for(i = 0; i < CHECK_ON_SIZE; i++){
          printf("%10.3f", jbi[1][i]);
        }

        if(compare(jbi[1], check_res, tab_size) == 0){
          printf("\nThe result with this method is not correct !\n");
        }

        printf("\n----------------------\n");
        printf("Total time :\t %13f ms\n", (double) accu * 1000.0);
        printf("Average time :\t %13f ms\n\n", 
          (double) (accu * 1000.0 / (nruns)));
      }
    }

    free(check_res);
    free(jbi[1]);
    free(jbi[0]);
    free(jbi);

}

/* --------*/
/*  Tests  */
/* --------*/

void djbi1d_sk_full_tiles_test(int n, int iters, double ** jbi, \
  struct benchscore * bsc){
  int steps = (iters / T_ITERS) + 1;
  int strips = (n/(T_WIDTH_DBL)) + 1;

  #ifdef DEBUG
    printf("iters : %i, n : %i -- %i steps, %i strips\n", steps, strips, 
      iters, n);
  #endif

  double * jbi_dashs = (double*) malloc(sizeof(double) * 
    ((strips + steps) * T_WIDTH_DBL));
  double * jbi_slashs = (double*) malloc(sizeof(double) * steps * 2 * T_ITERS);

  uint8_t ** tasks;

  if(jbi_dashs == NULL || jbi_slashs == NULL){
    fprintf(stderr, "Error while allocating 2D arrays for skewed_tiles\n");
    return;
  }

  for(int i = 0; i < n; i++){
    jbi_dashs[i] = jbi[0][i];
  }

  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  tasks = djbi1d_sk_full_tiles(strips, steps, jbi_dashs, jbi_slashs);
  clock_gettime( CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);

  #ifdef DEBUG
    if(task_index(tasks, strips, steps) > 0){
      printf("The task index for dependencies doesn't seem correct ...\n");
      for(int t = 0; t < steps; t ++){
        for(int i = 1; i < strips; i ++){
          printf("%2i", tasks[t][i]);
        }
        printf("\n");
      }
    }
  #endif
  int start_stripe_top = steps * T_WIDTH_DBL;
  for(int i = 0 ; i < n; i++){
    jbi[1][i] = jbi_dashs[i + start_stripe_top];
  }

  free(jbi_dashs);
  free(jbi_slashs);
  FREE_MX(tasks, steps)
}

void djbi1d_skewed_tiles_test(int n, int iters, double ** jbi, \
  struct benchscore * bsc){
  int steps = (iters / T_ITERS) + 1;
  int strips = (n/(T_WIDTH_DBL)) + 1;

  #ifdef DEBUG
    printf("iters : %i, n : %i -- %i steps, %i strips\n", steps, strips, 
      iters, n);
  #endif

  double * jbi_dashs = (double*) malloc(sizeof(double) * 
    ((strips + steps) * T_WIDTH_DBL));
  double * jbi_slashs = (double*) malloc(sizeof(double) * steps * 2 * T_ITERS);

  if(jbi_dashs == NULL || jbi_slashs == NULL){
    fprintf(stderr, "Error while allocating 2D arrays for skewed_tiles\n");
    return;
  }

  for(int i = 0; i < n; i++){
    jbi_dashs[i] = jbi[0][i];
  }

  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  djbi1d_skewed_tiles(strips, steps, jbi_dashs, jbi_slashs);
  clock_gettime( CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);

  int start_stripe_top = steps * T_WIDTH_DBL;
  for(int i = 0 ; i < n; i++){
    jbi[1][i] = jbi_dashs[i + start_stripe_top];
  }

  free(jbi_dashs);
  free(jbi_slashs);
}


/* Task functions :
  do_i0_t0 : starting task, bottom-left corner of the space-time matrix
  do_i0_t : left column, triangular tile
  do_i_t0 : first time iteration, bottom line
  do_i_t : regular task
  do_in_t : right column, triangular tile
*/

inline void do_i0_t0(double * dashs, double * slashs, int T, int I){

  ALLOC_LINES(l1, l2, (T_WIDTH_DBL + 1))

  uint8_t t,i;

  #ifdef DEBUG
    printf("First task T I, %i %i\n", T, I);
  #endif

  for(i = 1; i < T_WIDTH_DBL + 1; i++){
    l1[i] = dashs[i-1]; 
  }

  for(t = 0; t < T_ITERS * 2; t+=2){
    l1[0] = 0.0;
    int right = max(T_WIDTH_DBL - t/2, 1);
    for(i = 1; i < right; i++){
      l2[i] = (l1[i -1] + l1[i] + l1[i+1]) / 3.0 ;
    }
    slashs[t] = l1[right - 1];
    slashs[t + 1] = l1[right];
    memcpy(l1, l2,(T_WIDTH_DBL + 1)  * sizeof(double));
  }
  #ifdef DEBUG
    for(i = 0; i < CHECK_ON_SIZE; i++)
      printf("%10.3f", dashs[i]);
    printf("\n");
  #endif


  FREE_LINES(l1, l2)
}

inline void do_i0_t(double * dashs, double * slashs, int strpno, int Tt){

  ALLOC_LINES(l1, l2, (T_WIDTH_DBL + 1))

  uint8_t t,i;

  #ifdef DEBUG
    printf("Do do_i0_t %i %i\n", Tt, strpno - Tt);
  #endif

  for(i = 1; i < T_WIDTH_DBL; i++){
   l1[i] = dashs[strpno * T_WIDTH_DBL + i - 1];
  }

  for(int t= 0; t < 2 * T_ITERS; t += 2){
    l1[0] = 0.0;
    int right = max(T_WIDTH_DBL - t/2, 0);
    for(i = 1; i < right; i++){
      l2[i] = ((l1[i - 1] + l1[i] + l1[i + 1]) / 3.0);
    }
    slashs[Tt * 2 * T_ITERS + t] = l1[right - 1];
    slashs[Tt * 2 * T_ITERS + t + 1] = l1[right];

    memcpy(l1, l2,(T_WIDTH_DBL + 1)  * sizeof(double));
  }


  FREE_LINES(l1, l2)
}

inline void do_i_t(double * dashs, double * slashs, int strpno, int Tt){
  
  ALLOC_LINES(l1, l2, (T_WIDTH_DBL + 2))

  uint8_t t,i;

  #ifdef DEBUG
    if(Tt == ((DBG_ITER / T_ITERS ))){
      printf("Final line, task %i %i\n", Tt, strpno - Tt);
    }
  #endif

  // Load dash
  for(i = 2; i < T_WIDTH_DBL + 2; i++){
    l1[i] = dashs[strpno * T_WIDTH_DBL + i - 2];
  }

  for(t = 0; t < T_ITERS * 2; t+=2){
    // Load slash part
    l1[0] = slashs[Tt * 2 * T_ITERS + t];
    l1[1] = slashs[Tt * 2 * T_ITERS + t + 1];
    for(i = 2; i < T_WIDTH_DBL + 2; i++){
      l2[i] = ((l1[i -2] + l1[i - 1] + l1[i]) / 3.0) ;
    }
    // Write slash
    slashs[Tt * 2 * T_ITERS + t] = l1[T_WIDTH_DBL];
    slashs[Tt * 2 * T_ITERS + t + 1] = l1[T_WIDTH_DBL + 1];

    memcpy(l1, l2,(T_WIDTH_DBL + 2)  * sizeof(double));
  }
  // Write dash
  for(i = 2; i < T_WIDTH_DBL; i++){
    dashs[strpno * T_WIDTH_DBL + i - 2] = l2[i];
  }


  FREE_LINES(l1, l2)
}


inline void do_in_t(double * dashs, double * slashs, int strpno, int Tt){

  ALLOC_LINES(l1, l2, (T_WIDTH_DBL + 2))

  uint8_t t,i;

  #ifdef DEBUG
    if(Tt == ((DBG_SIZE / T_WIDTH_DBL ))){
      printf("Final task %i %i\n", Tt, strpno - Tt);
    } else {
      printf("Do do_in_t %i %i\n", Tt, strpno - Tt);
    }
  #endif

  for(t = 0; t < T_ITERS * 2; t+=2){
    // Load slash part
    l1[0] = slashs[Tt * 2 * T_ITERS + t];
    l1[1] = slashs[Tt * 2 * T_ITERS + t + 1];
    int r = (2 + t/2);
    l1[r] = 0.0;
    for(i = 2; i <= r; i++){
      l2[i] = ((l1[i -2] + l1[i - 1] + l1[i]) / 3.0) ;
    }
    memcpy(l1, l2,(T_WIDTH_DBL + 2)  * sizeof(double));
  }
  // Write dash
  for(i = 2; i < T_WIDTH_DBL; i++){
    dashs[strpno * T_WIDTH_DBL + i - 2] = l1[i];
  }

  FREE_LINES(l1, l2)
}