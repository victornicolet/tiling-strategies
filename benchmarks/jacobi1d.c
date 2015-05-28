#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include <inttypes.h>
#include <signal.h>
#include "jacobi1d.h"

static struct timespec tend;
static struct timespec tbegin;

static FILE * csv_file;

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

struct benchspec benchmarks[] = {
  {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap, check_tilable, 2, 1 << 13, 1 << 7},
  {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive, check_default, 2, 1 << 13, 1 << 7},
  {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test, check_tilable,
     2, 1 << 13, 1 << 7},
  {"JACOBI1D_SK_FULL_TILES", djbi1d_sk_full_tiles_test, check_tilable,
     2, 1 << 13, 1 << 7},
  {"JACOBI1D_SWAP_SEQ", djbi1d_swap_seq, check_default, 2, 1 << 13, 1 << 7},
  {"JACOBI1D_HALF_DIAMONDS", djbi1d_half_diamonds_test, check_low_iter,
   2, 1 << 18, 1 << 5}
};

void djbi1d_half_diamonds(int w, int iters, double * jbi){

  int base_w = 2 * iters ;
  int strips = (w / base_w);
  // Store the border between base-down pyramids and base-up pyramids
  int tmp_stride = 4 * iters - 2;
  double * tmp = (double *) malloc(sizeof(double) * tmp_stride * strips);

  int ti,t,i;
#ifdef DEBUG
  int prevr0 = -1;
  int * counters = (int *) calloc(w, sizeof(int));
  char ** viewtile = (char **) malloc(iters * sizeof(char*));
  for(i = 0; i < iters; i ++){
    viewtile[i] = (char *) malloc(w * sizeof(char));
    for(int ii = 0; ii < w; ii++){
      viewtile[i][ii] = 's';
    }
  }
  int track_cell = 5;
#endif


  // First loop : base down tiles
#ifndef SEQ
  #pragma omp parallel for schedule(static)
#endif
  for(ti = 0; ti < strips; ti ++){

    int tmp_pos = tmp_stride * ti;

    double li1[base_w], li0[base_w];
    // Initial values
    int l0 = max(ti * base_w, 0);
    int r0 = min(l0 + base_w, w);
    for(i = l0; i < r0; i ++){
      li0[i - l0] = jbi[i];
    }

    for(t = 1; t < iters; t ++){
      int l = max(l0 + t, 1);
      int r = min(l0 + base_w - t, w-1);
      // Fill the border-storing array
      tmp[tmp_pos + tmp_stride - 2*t]     = li0[r - l0];
      tmp[tmp_pos + tmp_stride - 2*t - 1] = li0[r - l0 - 1];
      tmp[tmp_pos + 2*(t-1)]              = li0[l - l0];
      tmp[tmp_pos + 2*(t-1) + 1]          = li0[l - l0 + 1];

      for(i = l; i < r; i ++){
        li1[i - l0] = (li0[i-1-l0] + li0[i - l0] + li0[i+1-l0]) / 3.0;
        #ifdef DEBUG
          viewtile[t-1][i] = 'X';
          counters[i]++;
          if(i == track_cell){
            printf("%i, %i : %10.3f\n", i, t, li0[i]);
          }
        #endif
      }
      for(i = l; i < r; i ++){
        li0[i - l0] = li1[i - l0];
      }
    }
  }
#ifdef DEBUG
  for(int tt = iters - 1; tt >= 0; tt --){
    for(int ii = 0; ii < 80; ii++){
      printf("%c", viewtile[tt][ii]);
    }
    printf("\n");
  }
#endif
 // --------------------------------------
#ifdef GDB_DEBUG
  raise(SIGABRT);
#endif
  // Second loop : tip down tiles
#ifndef SEQ
  #pragma omp parallel for schedule(static)
#endif
  for(ti = 0; ti < strips + 1; ti ++){

    int tmp_pos = tmp_stride * ti;
    int x0 = ti * base_w;
    int l0 = max(x0 - iters - 1, 0);
    int r0 = min(x0 + iters, w);
    double li1[base_w + 2], li0[base_w + 2];

    for(t = 0; t < iters; t ++){
      int l = max(x0 - (t+1), 1);
      int r = min(x0 + t + 1, w);
      // Load from the border-storing array
      li0[r - l0 - 1]     =
        tmp[min(tmp_pos + 2 * (t - 1), tmp_stride * strips)];
      li0[r - l0] =
        tmp[min(tmp_pos + 2 * (t - 1) - 1, tmp_stride * strips)];
      li0[l - l0 - 1] =
        tmp[max(tmp_pos - 2 * (t - 1), 0)];
      li0[l - l0] =
        tmp[max(tmp_pos - 2 * (t - 1) + 1, 0)];

      for(i = l; i < r; i ++){
        li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
        #ifdef DEBUG
          if(viewtile[t][i] == 's'){
            viewtile[t][i] = '.';
          } else {
            viewtile[t][i] = 'o';
          }
          counters[i]++;
          if(i == track_cell){
            printf("%i, %i : %10.3f\n", i, t, li0[i]);
          }
        #endif
      }
      for(int ii = 0; ii < base_w + 2; ii ++){
        li0[ii] = li1[ii];
      }
    }
  // --------------------------------------
#ifdef GDB_DEBUG
  raise(SIGABRT);
#endif

  // Debug general tile layout
  #ifdef DEBUG
    #ifdef SEQ
      if(prevr0 > l0 + 1 && prevr0 != -1){
        fprintf(stderr, "prevr0 : %i r0 : %i l0 : %i \n", prevr0, r0, l0);
        fprintf(stderr, "Error ! Top of tile %i overlaps with previous tile ! \n \
          Aborting...\n",
         ti);
        return;
      } else if((prevr0 - l0 - 1) > 1 && prevr0 != -1){
        fprintf(stderr, "prevr0 : %i r0 : %i l0 : %i \n", prevr0, r0, l0);
        fprintf(stderr, "Error ! Top of tile %i too far from previous tile ! \n \
          Aborting...\n",
         ti);
        return;
      }
      prevr0 = r0;
    #endif
  #endif

    // Copy back to memory
    for(i = l0 + 1; i < r0; i++){
      jbi[i] = li0[i - l0];
    }
  }
  // Check counters
  #ifdef DEBUG
    int pv = 0, counting = 0, prints = 0;
    for(i = 1; i < w - 1; i ++){
      if(counters[i] != iters){
        if(pv != i - 1){
          fprintf(stderr, "Error : bad operations count from cell %i ", i);
          counting = 1;
        }
        pv = i;
      }
      if(i != pv && counting == 1){
        fprintf(stderr, "to cell %i. ", i );
        fprintf(stderr, "[%i operations]\n", counters[pv]);
        counting = 0;
        prints ++;
      }
      if(prints > 3){
        fprintf(stderr, "Too much cells with bad operation counts..\n");
        break;
      }
    }
    printf("\n");
    for(int tt = iters - 1; tt >= 0; tt --){
      for(int ii = 0; ii < 80; ii++){
        printf("%c", viewtile[tt][ii]);
      }
      printf("\n");
    }
  #endif

}

  /* In this algorithm we work only on full parallelograms : no partial tiles.
  We are interested in performance measures, and understanding, rather
  than corectness here */
uint8_t ** djbi1d_sk_full_tiles(int num_strips, int num_steps, double * dashs, \
  double * slashs){

  int i,t;

  uint8_t ** taskdep_index = malloc(num_steps * sizeof ** taskdep_index);
  for(i = 0; i < num_strips; i++){
    taskdep_index[i] = malloc(num_steps * sizeof * taskdep_index[i]);
  }


#ifndef SEQ
  #pragma omp parallel
#endif
  {
#ifndef SEQ
    #pragma omp single
#endif
    {
      for(t = 0; t < num_steps; t ++){
        for(i = 1; i < num_strips; i++){
          taskdep_index[t][i] = 0;
          int strpno = i + t;
          // Stratup task
          if(t == 0 && i == 1){
#ifndef SEQ
            #pragma omp task firstprivate(i,t) \
            depend(out : taskdep_index[t][i])
#endif
            {
              do_i_t(dashs, slashs, 1, 0);
              taskdep_index[t][i] ^= 1;
            }
          }else if(t == 0 && i > 1){
#ifndef SEQ
            #pragma omp task firstprivate(i,t) \
            depend(out : taskdep_index[t][i])\
            depend(in : taskdep_index[i-1][0])
#endif
            {
              do_i_t(dashs, slashs, strpno, 0);
              taskdep_index[t][i] ^= 1;
            }
          } else {
#ifndef SEQ
            #pragma omp task firstprivate(i,t) \
            depend(out : taskdep_index[t][i]) \
            depend(in  : taskdep_index[t][i-1], taskdep_index[t-1][i+1])
 #endif
            {
              do_i_t(dashs, slashs, strpno, t);
              taskdep_index[t][i] ^= 1;
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

    #ifdef DEBUG_PARALLEL
      int * tsk = (int*) malloc(sizeof(int) * (strips + 1) * steps);
    #endif
    int Ti, Tt;
#ifndef SEQ
  #pragma omp parallel
  #pragma omp master
#endif
    {
      for(Tt = 0; Tt < steps; Tt++){
        for(Ti = 0; Ti < strips + 1; Ti++){
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
              #ifdef DEBUG_PARALLEL
                tsk[0] = 1;
              #endif
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
              #ifdef DEBUG_PARALLEL
                if(tsk[Ti - 1] != 1 ){
                  printf("Unsatisified dependency !\n");
                }
                tsk[Tt * steps + Ti] = 1;
              #endif
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
              #ifdef DEBUG_PARALLEL
                if(tsk[(Tt - 1) * steps + Ti] != 1 ){
                  printf("Unsatisified dependency !\n");
                }
                tsk[Tt * steps + Ti] = 1;
              #endif
              do_i0_t(dashs, slashs, strpno, Tt);
            }
          } else if(Ti == strips){
#ifndef SEQ
  #pragma omp task firstprivate(Ti,Tt) \
    depend(in: tasks[Ti-1][Tt]) \
    depend(out : tasks[Ti][Tt])
#endif
            {
              #ifdef DEBUG_PARALLEL
                if(tsk[Tt * steps + Ti - 1] != 1 ){
                  printf("Unsatisified dependency !\n");
                }
                tsk[Tt * steps + Ti] = 1;
              #endif
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
              #ifdef DEBUG_PARALLEL
                if(tsk[Tt * steps + Ti - 1] != 1 || \
                  tsk[(Tt-1) * steps + (Ti+1)] != 1){
                  printf("Unsatisified dependency !\n");
                }
                tsk[Tt * steps + Ti] = 1;
              #endif
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

  for( Tt = 0; Tt <= jbi_iters/T_ITERS; Tt ++){
#ifndef SEQ
    #pragma omp parallel for schedule(static)
#endif
    for(Ti = 0; Ti <= (n / T_WIDTH_DBL_OVERLAP); Ti ++){

      double * lvl1 = (double *) malloc(sizeof(double) * (T_WIDTH_DBL_OVERLAP+
          T_ITERS * 2));
      double * lvl0 = (double *) malloc(sizeof(double) * (T_WIDTH_DBL_OVERLAP+
          T_ITERS * 2));

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
    l2[0] = l1[0];
    l2[n - 1] = l1[n - 1];
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

void jbi_init(double **jbi, int n) {
  int j;

  for (j = 0; j < n; j++) {
    jbi[0][j] = fabs(cos((double) j )) * (1 << 8);
    jbi[1][j] = 0;
  }
  jbi[0][0] = 0.0;
  jbi[0][n-1] = 0.0;
}

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
    int tab_size, iterations;
    if(argc == 5){
      tab_size = 1 << atoi(argv[3]);
      iterations = 1 << atoi(argv[4]);
    } else {
      #ifdef DEBUG
        tab_size = DBG_SIZE;
        iterations = DBG_ITER;
      #endif
    }

    #ifndef DEBUG
      tab_size = ((tab_size - 1) / (2*T_ITERS - 1)) * (2*T_ITERS - 1);
    #endif

    char *benchmask = argv[2];
    if(strlen(benchmask) != nbench){
      printf("Error : not a valid mask ! Your mask must be %i bits long\n",
        nbench);
      return -1;
    }

    #ifdef DEBUG
      fprintf(csv_file, "Where ?; Data\n");
      fprintf(csv_file, "Input;\n");
    #endif

    printf("\n");
    double accu;

    for(int bs = 0; bs < nbench; bs++){
      if (benchmask[bs] == '1') {

        struct benchscore score[nruns + 1];
        accu = 0.0;
        if(argc < 5){
          #ifndef DEBUG
            tab_size = benchmarks[bs].size;
            iterations = benchmarks[bs].iters;
          #endif
        }

        if(benchmarks[bs].checkfunc(tab_size, iterations) < 0){
          fprintf(stderr, "Error : argument incompatible with variant\n");
          fprintf(stderr, "Iterations : %i \t Width : %i\n",
            iterations, tab_size);
          fprintf(stderr, "Variant : %s\n", benchmarks[bs].name);
          #ifndef DEBUG
            continue;
          #endif
        }

        double ** jbi = (double **) malloc(sizeof(double) * 2);
        for(i = 0; i < 2; i++){
          jbi[i] = (double *) malloc(sizeof(double) * tab_size);
        }
        // Get the correct result
        jbi_init(jbi, tab_size);
        double * check_res = (double *) malloc(sizeof(double) * tab_size);
        struct benchscore bsc;
        djbi1d_swap_seq(tab_size, iterations, jbi, &bsc);
        for(i = 0; i < tab_size; i++) check_res[i] = jbi[1][i];

        jbi_init(jbi, tab_size);
        printf("Input : \n");
        for(i = 0; i < DISPLAY_SIZE; i++){
          printf("%10.3f", jbi[0][i]);
        }
        printf("\n");
        #ifdef DEBUG
          for(i = 0; i < 1 << 8; i++){
            fprintf(csv_file, "%10.3f ;", jbi[0][i]);
          }
          fprintf(csv_file, "\n");
        #endif

        for(iter = 0; iter < nruns + 1; iter++){
          jbi_init(jbi, tab_size);
          score[iter].name = benchmarks[bs].name;
          benchmarks[bs].variant(tab_size, iterations, jbi, &score[iter]);

          if(iter > 0) {
            printf("%s : Run %i ...", score[iter].name, iter );
            printf("\t\t %13f ms\n", score[iter].wallclock * 1000.0 );
            accu += score[iter].wallclock ;
          }
        }
        printf("\n------------- %s ---------\n", benchmarks[bs].name);
        printf("Result: \n");
        for(i = 0; i < DISPLAY_SIZE; i++){
          printf("%10.3f", jbi[1][i]);
        }

        if(compare(jbi[1], check_res, tab_size) == 0){
          printf("\nThe result with this method is not correct ! ");
          printf("This should be the correct result :\n");
          for(i = 0; i < DISPLAY_SIZE; i++){
            printf("%10.3f", check_res[i]);
          }
        }
        printf("\n----------------------\n");
        printf("Total time :\t %13f ms\n", (double) accu * 1000.0);
        printf("Average time :\t %13f ms\n\n",
          (double) (accu * 1000.0 / (nruns)));

        free(check_res);
        free(jbi[1]);
        free(jbi[0]);
        free(jbi);
      }
    }
    #ifdef SEQ
          printf("WARNING : SEQ defined\n");
    #endif
    #ifdef DEBUG_PARALLEL
          printf("WARNING : DEBUG_PARALLEL defined\n");
    #endif
}

/* --------*/
/*  Tests  */
/* --------*/

void djbi1d_half_diamonds_test(int n, int iters, double ** jbi, \
  struct benchscore * bsc){
  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  djbi1d_half_diamonds(n, iters, jbi[0]);
  clock_gettime( CLOCK_MONOTONIC, &tend);
  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
  memcpy(jbi[1], jbi[0], sizeof(double) * n);
}


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
  do_i0_t0 : starting task, botile_tom-left corner of the space-time matrix
  do_i0_t : left column, triangular tile
  do_i_t0 : first time iteration, botile_tom line
  do_i_t : regular task
  do_in_t : right column, triangular tile
*/

inline void
do_i0_t0(double * dashs, double * slashs, int tile_t, int tile_i)
{

  double * l1 = alloc_line(T_WIDTH_DBL + 1);
  double * l2 = alloc_line(T_WIDTH_DBL + 1);

  uint8_t t,i;

  #ifdef DEBUG
    printf("First task T I, %i %i\n", tile_t, tile_i);
  #endif

  for (i = 1; i < T_WIDTH_DBL + 1; i++){
    l1[i] = dashs[i-1];
  }

  for (t = 0; t < T_ITERS * 2; t+=2){
    l1[0] = 0.0;
    int right = max(T_WIDTH_DBL - t/2, 1);

    memcpy(l2, l1,(T_WIDTH_DBL + 1)  * sizeof(double));
    for (i = 1; i < right; i++){
      l2[i] = (l1[i -1] + l1[i] + l1[i+1]) / 3.0 ;
    }

    slashs[t] = l1[right - 1];
    slashs[t + 1] = l1[right];
    memcpy(l1, l2,(T_WIDTH_DBL + 1)  * sizeof(double));
  }
  #ifdef DEBUG
    for(i = 0; i < DISPLAY_SIZE; i++)
      printf("%10.3f", dashs[i]);
    printf("\n");
  #endif


  free(l1);
  free(l2);
}

inline void
do_i0_t(double * dashs, double * slashs, int strpno, int tile_t)
{

  double * l1 = alloc_line(T_WIDTH_DBL + 1);
  double * l2 = alloc_line(T_WIDTH_DBL + 1);

  uint8_t t,i;

  #ifdef DEBUG
    printf("Do do_i0_t %i %i\n", tile_t, strpno - tile_t);
  #endif

  for(i = 1; i < T_WIDTH_DBL + 1; i++){
   l1[i] = dashs[strpno * T_WIDTH_DBL + i - 1];
  }

  for(t = 0; t < 2 * T_ITERS; t += 2){
    l1[0] = 0.0;
    int right = max(T_WIDTH_DBL - t/2, 0);

    memcpy(l2, l1,(T_WIDTH_DBL + 1)  * sizeof(double));
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

    memcpy(l2, l1,(T_WIDTH_DBL + 2)  * sizeof(double));
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


  free(l1);
  free(l2);
}


inline void do_in_t(double * dashs, double * slashs, int strpno, int tile_t){

  double * l1 = alloc_line(T_WIDTH_DBL + 2);
  double * l2 = alloc_line(T_WIDTH_DBL + 2);

  uint8_t t,i;

  #ifdef DEBUG
    if(tile_t == ((DBG_SIZE / T_WIDTH_DBL ))){
      printf("Final task %i %i\n", tile_t, strpno - tile_t);
    } else {
      printf("Do do_in_t %i %i\n", tile_t, strpno - tile_t);
    }
  #endif

  for(t = 0; t < T_ITERS * 2; t+=2){
    // Load slash part
    l1[0] = slashs[tile_t * 2 * T_ITERS + t];
    l1[1] = slashs[tile_t * 2 * T_ITERS + t + 1];
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

  free(l1);
  free(l2);
}
