#define _POSIX_C_SOURCE 200112L

#include <immintrin.h>
#include <inttypes.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "jacobi1d.h"

static struct timespec tend;
static struct timespec tbegin;

/* Functions describing different tasks in the computation
  Always inlined in the main body */
void do_i_t(double *, double *, int , int)
    __attribute__((always_inline));
void do_i0_t(double *, double *, int, int)
    __attribute__((always_inline));
void do_i0_t0(double *, double *, int, int)
    __attribute__((always_inline));
void do_in_t(double *, double *, int, int)
    __attribute__((always_inline));

struct benchspec benchmarks[] = {
  {"JACOBI1D_OMP_OVERLAP", djbi1d_omp_overlap, check_tilable,
        2, 1 << 13, 1 << 8},
  {"JACOBI1D_OMP_NAIVE", djbi1d_omp_naive, check_default,
        2, 1 << 13, 1 << 8},
  {"JACOBI1D_SKEWED_TILES", djbi1d_skewed_tiles_test, check_tilable,
        2, 1 << 13, 1 << 5},
  {"JACOBI1D_SK_FULL_TILES", djbi1d_sk_full_tiles_test, check_tilable,
        2, 1 << 13, 1 << 5},
  {"JACOBI1D_SWAP_SEQ", djbi1d_sequential, check_default,
        2, 1 << 13, 1 << 5},
  {"JACOBI1D_HALF_DIAMONDS", djbi1d_half_diamonds_test, check_low_iter,
        2, 1 << 13, 1 << 5}
};

/*
* task[i][j] is set to 1 if the task on time step i and column j has been
* executed. This function checks if there is no task that has been executed
* without its dependencies being statisfied.
*/
int
task_index(uint8_t ** tasks, int num_strips, int num_steps)
{
  int i,t;

  for (i = 1; i < num_strips; i++) {
    if (tasks[0][i] == 0 && tasks[0][i + 1] == 1) return -1;
  }

  for (t = 1; t < num_steps; t++) {
    for (i = 1; i < num_strips; i ++) {
      if ((tasks[i + 1][t - 1] == 0 || tasks[i - 1][t] == 0) &&
         tasks[i][t] == 1)
            return -1;
    }
  }
  return 1;
}

int
check_low_iter(int pb_size, int num_stencil_iters)
{
  if ((num_stencil_iters >= 16) && (num_stencil_iters <= 64) &&
    (pb_size > num_stencil_iters *(1 << 3))) {
      return 1;
  } else {
      return -1;
  }
}

int
check_tilable(int pb_size, int num_stencil_iters)
{
  if (pb_size > (T_WIDTH_DBL << 2) &&
    num_stencil_iters > 2*T_ITERS) {
    return 1;
  } else {
    return -1;
  }
}

int
check_default(int pb_size, int num_stencil_iters)
{
  if ( pb_size > num_stencil_iters && num_stencil_iters > 2) {
    return 1;
  } else {
    return -1;
  }
}

/* Different versions of jacobi1d have been implemented here.
*   - with half - diamonds : useful if we have few iterations, it allows
*      parallelism an a concurrent start.
*   - skewed tiles, with tasks
*   - with overlapping tiles
*   - sequential with swapping
*/

void djbi1d_half_diamonds(int pb_size, int num_iters, double * jbi) {

  int tile_no, t, i, ii, tmp_pos;
  // Tile bounds
  int l, l0, r, r0, x0;

  int tile_base_sz = 2 * num_iters ;
  int num_tiles = (pb_size / tile_base_sz);
  // Store the border between base-down pyramids and base-up pyramids
  int tmp_stride = 4 * num_iters - 2;
  double * tmp = malloc(tmp_stride * num_tiles * sizeof(*tmp));

#ifdef DEBUG
  int * counters = calloc(pb_size, sizeof(*counters));
  char ** viewtile = malloc(num_iters * sizeof(*viewtile));

  for (i = 0; i < num_iters; i ++) {
    viewtile[i] = malloc(pb_size * sizeof(*viewtile[i]));
    for (ii = 0; ii < pb_size; ii ++) {
      viewtile[i][ii] = 's';
    }
  }
  int track_cell = 2;
#endif


/* First loop : base down tiles */
#ifndef SEQ
  #pragma omp parallel for schedule(static) shared(tmp) \
   private(tmp_pos, l0, r0, l, r, x0)
#endif
  for (tile_no = 0; tile_no < num_tiles; tile_no ++) {

    tmp_pos = tmp_stride * tile_no;

    double * li1 = alloc_line(tile_base_sz + 2);
    double * li0 = alloc_line(tile_base_sz + 2);

    /* Initial values */
    l0 = max(tile_no * tile_base_sz, 0);
    r0 = min(l0 + tile_base_sz, pb_size - 1);
    for (i = l0; i <= r0; i ++) {
      li0[i - l0] = jbi[i];
      li1[i - l0] = 0.0;
    }
    for(i = r0 + 1; i < tile_base_sz + l0; i++){
      li0[i - l0] = 0.0;
      li1[i - l0] = 0.0;
    }

    for (t = 1; t < num_iters; t ++) {
      l = max(l0 + t, 1);
      r = min(l0 + tile_base_sz - t, pb_size);
/* The border of the pyramids needs to be stored for further computation */
      tmp[tmp_pos + tmp_stride - 2 * t - 1]     = li0[r - l0];
      tmp[tmp_pos + tmp_stride - 2 * t - 2]     = li0[r - l0 - 1];
      tmp[tmp_pos + 2*(t-1)]                    = li0[l - l0 - 1];
      tmp[tmp_pos + 2*(t-1) + 1]                = li0[l - l0];

      for (i = l; i < r; i ++) {
        li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
#ifdef DEBUG
          if((i == r - 1) || (i == r - 2) || (i == l) || (i == l + 1)){
            viewtile[t-1][i] = 'x';
          } else {
            viewtile[t-1][i] = 'X';
          }
          counters[i] ++;
          if (i == track_cell) {
            printf("%i, %i : %10.3f\n", i, t, li0[i]);
          }
#endif
      }
      for (i = l; i < r; i ++) {
        li0[i - l0] = li1[i - l0];
      }
    }
  }

/* Second loop : tip down tiles */
#ifndef SEQ
  #pragma omp parallel for schedule(static) shared(tmp) \
  private(tmp_pos, l0, r0, l, r, x0)
#endif
  for (tile_no = 0; tile_no < num_tiles + 1; tile_no ++) {

    tmp_pos = tmp_stride * tile_no;
    x0 = tile_no * tile_base_sz;
    l0 = max(x0 - num_iters - 1, 0);
    r0 = min(x0 + num_iters, pb_size);

    double * li1 = alloc_line(tile_base_sz + 2);
    double * li0 = alloc_line(tile_base_sz + 2);

    for (t = 0; t < num_iters; t ++) {
      l = max(x0 - (t+1), 1);
      r = min(x0 + t, pb_size);
      /* Load from the border-storing array */
       li0[max(r - l0 - 1, 0)]  =
        tmp[max(min(tmp_pos + 2 * t, tmp_stride * num_tiles - 1), 0)];
       li0[r - l0] =
        tmp[max(min(tmp_pos + 2 * t + 1, tmp_stride * num_tiles - 1), 0)];
       li0[l - l0 - 1] =
        tmp[min(max(tmp_pos - 2 * t - 1, 0), tmp_stride * num_tiles - 1)];
       li0[l - l0] =
        tmp[min(max(tmp_pos - 2 * t - 2, 0), tmp_stride * num_tiles - 1)];


      for (i = l; i <= r; i ++) {
        li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
#ifdef DEBUG
          if (viewtile[t][i] == 's') {
            viewtile[t][i] = '.';
          } else {
            viewtile[t][i] = 'o';
          }
          counters[i]++;
          if (i == track_cell) {
            printf("%i, %i : %10.3f\n", i, t, li0[i - l0]);
          }
#endif
      }

      for (ii = 0; ii < tile_base_sz + 2; ii ++) {
        li0[ii] = li1[ii];
      }
    }
    /* Copy back to memory */
    for (i = l0 + 1; i < r0; i++) {
      jbi[i] = li0[i - l0];
    }
    free(li0);
    free(li1);


  }
  free(tmp);

/* ---------------------------*/
/* Only debug below this line */

#ifdef DEBUG
  int pv = 0, counting = 0, prints = 0;
  for (i = 1; i < pb_size - 1; i ++) {
    if (counters[i] != num_iters) {
      if (pv != i - 1) {
        fprintf(stderr, "Error : bad operations count from cell %i ", i);
        counting = 1;
      }
      pv = i;
    }
    if (i != pv && counting == 1) {
      fprintf(stderr, "to cell %i. ", i );
      fprintf(stderr, "[%i operations]\n", counters[pv]);
      counting = 0;
      prints ++;
    }
    if (prints > 3) {
      fprintf(stderr, "Too much cells with bad operation counts..\n");
      break;
    }
  }
  printf("\n");
  for (int tile_t = num_iters - 1; tile_t >= 0; tile_t --) {
    printf("%i \t- ", tile_t);
    for (int ii = 0; ii < 80; ii++) {
      if(viewtile[tile_t][ii] == 'x'){
        printf("%s%c%s", KBLU, viewtile[tile_t][ii], KRESET);
      } else {
        printf("%c", viewtile[tile_t][ii]);
      }
    }
    free(viewtile[tile_t]);
    printf("\n");
  }
  free(viewtile);
  free(counters);
#endif

}



/* In this algorithm we work only on full parallelograms : no partial tiles.
 * We are interested in performance measures, and understanding, rather
 * than corectness here
 */
uint8_t **
djbi1d_sk_full_tiles(int num_strips, int num_steps, double * dashs,
  double * slashs)
{

  int i,t;

  uint8_t ** taskdep_index = malloc(num_steps * sizeof(*taskdep_index));
  for (i = 0; i < num_strips; i ++) {
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
      for (t = 0; t < num_steps; t ++) {
        for (i = 1; i < num_strips; i++) {
          taskdep_index[t][i] = 0;
          int strpno = i + t;
          if (t == 0 && i == 1) {
#ifndef SEQ
            #pragma omp task firstprivate(i,t)                                \
            depend(out : taskdep_index[t][i])
#endif
            {
              do_i_t(dashs, slashs, 1, 0);
              taskdep_index[t][i] ^= 1;
            }
          }else if (t == 0 && i > 1) {
#ifndef SEQ
            #pragma omp task firstprivate(i,t)                                \
            depend(out : taskdep_index[t][i])                                 \
            depend(in : taskdep_index[i-1][0])
#endif
            {
              do_i_t(dashs, slashs, strpno, 0);
              taskdep_index[t][i] ^= 1;
            }
          } else {
#ifndef SEQ
            #pragma omp task firstprivate(i,t)                                \
            depend(out : taskdep_index[t][i])                                 \
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

  return taskdep_index;
}

void
djbi1d_skewed_tiles(int num_strips, int num_steps, double * dashs,
  double * slashs)
{
/* In this version we assume T_ITERS = T_WIDTH_DBL so we have to make a
 *   difference between regular tiles ( parallelogram-shaped ones) and
 *   triangular tiles, but the patile_tern is quite regular and dependencies are
 *   straightforward at the boundaries of the domain.
 *   ---------------
 *
 *  num_steps = (iters / T_ITERS) + 1;
 *   num_strips = (n / T_WIDTH_DBL) + 1 ;
 */

    /* Unused except in omp task pragmas */
    uint8_t taskdep_index[num_strips+1][num_steps] __attribute__((unused));

    #ifdef DEBUG_PARALLEL
      int * tsk = (int*) malloc(sizeof(int) * (num_strips + 1) * num_steps);
    #endif
    int tile_i, tile_t;
#ifndef SEQ
  #pragma omp parallel
  #pragma omp master
#endif
    {
      for (tile_t = 0; tile_t < num_steps; tile_t++) {
        for (tile_i = 0; tile_i < num_strips + 1; tile_i++) {
          /* Strip number */
          int strpno = (tile_i + tile_t);

          if ( tile_t == 0 && tile_i == 0) {
#ifndef SEQ
            #pragma omp task firstprivate(tile_i,tile_t)                      \
              depend(out : taskdep_index[0][0])
#endif
            {
              #ifdef DEBUG_PARALLEL
                tsk[0] = 1;
              #endif
              do_i0_t0(dashs, slashs, tile_i, tile_t);
            }
          } else if (tile_t == 0 && tile_i < num_strips - 1) {
            /* Bottom tiles : only left-to-right dependencies + top out */
#ifndef SEQ
            #pragma omp task firstprivate(tile_i,tile_t)                      \
              depend(in : taskdep_index[tile_i-1][0])                         \
              depend(out : taskdep_index[tile_i][tile_t])
#endif
            {
              #ifdef DEBUG_PARALLEL
                if (tsk[tile_i - 1] != 1 ) {
                  printf("Unsatisified dependency !\n");
                }
                tsk[tile_t * num_steps + tile_i] = 1;
              #endif
              do_i_t(dashs, slashs, strpno, 0);
            }

          } else if (tile_i == 0 && tile_t > 0) {
            /* Left edge tile : triangular tiles
             * Only one in dependency, one out
             * ( here we assume T_ITERS == T_WIDTH_DBL )
             */
#ifndef SEQ
            #pragma omp task firstprivate(tile_i,tile_t)                      \
              depend(in : taskdep_index[1][tile_t-1])                         \
              depend(out : taskdep_index[tile_i][tile_t])
#endif
            {
              #ifdef DEBUG_PARALLEL
                if (tsk[(tile_t - 1) * num_steps + tile_i] != 1 ) {
                  printf("Unsatisified dependency !\n");
                }
                tsk[tile_t * num_steps + tile_i] = 1;
              #endif
              do_i0_t(dashs, slashs, strpno, tile_t);
            }
          } else if (tile_i == num_strips) {
#ifndef SEQ
            #pragma omp task firstprivate(tile_i,tile_t)                      \
              depend(in: taskdep_index[tile_i-1][tile_t])                     \
              depend(out : taskdep_index[tile_i][tile_t])
#endif
            {
              #ifdef DEBUG_PARALLEL
                if (tsk[tile_t * num_steps + tile_i - 1] != 1 ) {
                  printf("Unsatisified dependency !\n");
                }
                tsk[tile_t * num_steps + tile_i] = 1;
              #endif
              do_in_t(dashs, slashs, strpno, tile_t);
            }

          }else{
            /* Regular tile two in and out dependencies */
#ifndef SEQ
            #pragma omp task firstprivate(tile_i,tile_t)                      \
              depend(in : taskdep_index[tile_i-1][tile_t],                    \
                        taskdep_index[tile_i+1][tile_t-1])                    \
              depend(out : taskdep_index[tile_i][tile_t])
#endif
            {
              #ifdef DEBUG_PARALLEL
                if (tsk[tile_t * num_steps + tile_i - 1] != 1 || \
                  tsk[(tile_t-1) * num_steps + (tile_i+1)] != 1) {
                  printf("Unsatisified dependency !\n");
                }
                tsk[tile_t * num_steps + tile_i] = 1;
              #endif
              do_i_t(dashs, slashs, strpno, tile_t);
            }
          }
        }
      }
  }
}



void
djbi1d_diamond_tiles(int n,int num_stencil_iters, double ** jbi,
  struct benchscore * bsc)
{
  int r, l, bot, top;

  int stg = (num_stencil_iters / T_ITERS_DIAM) + 1;
  int strp = (n / T_WIDTH_DBL_DIAM);
  for (long tile_i = -1; tile_i < strp + 1; tile_i++) {
    for (long tile_t = 0; tile_t < stg; tile_t++) {
      /* tile_ile height */
      bot = max(tile_t * T_ITERS_DIAM, 1);
      top = min((tile_t + 1) * T_ITERS_DIAM, num_stencil_iters);
      for (int t = bot; t < top; t++) {
        /* Line boundaries */
        l = max(tile_i * T_WIDTH_DBL_DIAM - (t - bot) , 1);
        r = min((tile_i + 1) * T_WIDTH_DBL_DIAM - (t - bot), n - 1);
        for (int i = l; i < r; i++) {
          JBI1D_STENCIL_T(jbi);
        }
      }
    }
  }
}

void
djbi1d_omp_naive(struct args_dimt args, double * jbi_in, double * jbi_out,
  struct benchscore * bsc)
{
  int n = args.width;
  int num_stencil_iters = args.iters;
  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  int t,i;
  double * l1 = (double *) aligned_alloc(CACHE_LINE_SIZE, sizeof(double) * n);
  double * l2 = (double *) aligned_alloc(CACHE_LINE_SIZE, sizeof(double) * n);
  memcpy(l1, jbi_in, n * sizeof(double));

  clock_gettime( CLOCK_MONOTONIC, &tbegin);

  for (t = 0; t < num_stencil_iters; t++) {
#ifdef SEQ
   #pragma omp parallel for schedule(static)
#endif
    for (i = 1; i < n - 1; i++) {
      JBI1D_STENCIL(l2, l1);
    }
    memcpy(l1, l2, n * sizeof(double));
  }

  clock_gettime( CLOCK_MONOTONIC, &tend);

  for (int i = 0; i < n; i++) {
    jbi_out[i] = l1[i];
  }

  free(l1);
  free(l2);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
}

void
djbi1d_omp_overlap(struct args_dimt args, double * jbi_in, double * jbi_out,
  struct benchscore * bsc)
{
  int pb_size = args.width, num_stencil_iters = args.iters;
  int tile_i, tile_t, t, i;
  int tile_base_sz = T_WIDTH_DBL_OVERLAP + T_ITERS * 2;

  clock_gettime(CLOCK_MONOTONIC, &tbegin);

  for ( tile_t = 0; tile_t <= num_stencil_iters/T_ITERS; tile_t ++) {
#ifndef SEQ
    #pragma omp parallel for schedule(static)
#endif
    for (tile_i = 0; tile_i <= (pb_size / T_WIDTH_DBL_OVERLAP); tile_i ++) {

      double * lvl1 = malloc(tile_base_sz * sizeof(*lvl1));
      double * lvl0 = malloc(tile_base_sz * sizeof(*lvl0));


      /* Compute tile bounds */
      int bot = max((T_ITERS * tile_t), 1);
      int top = min((T_ITERS * (tile_t + 1)), num_stencil_iters);
      int h = top - bot;
      int l0 = max((T_WIDTH_DBL_OVERLAP * tile_i), 0);
      int r0 = min((T_WIDTH_DBL_OVERLAP * (tile_i + 1)), pb_size);
      int l = max((l0 - h ), 0);
      int r = min((r0 + h), pb_size);
      {
        /* Read tile base */
        for (i = l ; i < r ; i++ ) {
          lvl0[i- l] = jbi_in[i];
          lvl1[i- l] = 0.0;
        }

        for (t = 0 ; t < h; t++) {

          int lt = max(l0 - t , 1);
          int rt = min(r0 + t, pb_size);

          for (i = lt ; i < rt; i++) {
            lvl1[i - l] = (lvl0[i - l - 1] + lvl0[i - l] + lvl0[i - l + 1])
              / 3.0;
          }
          memcpy(lvl0, lvl1,
            (T_WIDTH_DBL_OVERLAP + T_ITERS * 2) * sizeof(double));

#ifdef DEBUG
          if (tile_i == 0) {
              fprintf(csv_file, "Left column ; %i ;%i", tile_t, t);
              for (i = 0; i < 8 ; i++) {
                fprintf(csv_file, ";%10.3f", lvl1[i] - jbi_in[i + l]);
              }
              fprintf(csv_file, "\n");
          }
#endif
        }

        /* Write tile top */
        for (i = l0 ; i < r0 ; i++ ) {
          jbi_out[i] = lvl0[i-l];
        }
      }
      free(lvl1);
      free(lvl0);
    }
/* Implicit barrier here, when all chunks of the loops are finished,
 * we copy the resulting data in the "source"
 */
    memcpy(jbi_in, jbi_out, pb_size * sizeof(double));
  }

  clock_gettime( CLOCK_MONOTONIC, &tend);



  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
}


void
djbi1d_sequential(struct args_dimt args, double * jbi_in, double * jbi_out,
 struct benchscore * bsc)
{
  int num_stencil_iters = args.iters, n = args.width;
  /* Boundaries initial condition */
  int t,i;
  double * l1 = (double *) aligned_alloc(CACHE_LINE_SIZE, sizeof(*l1) * n);
  double * l2 = (double *) aligned_alloc(CACHE_LINE_SIZE, sizeof(*l2) * n);
  memcpy(l1, jbi_in, n * sizeof(*jbi_in));

  clock_gettime( CLOCK_MONOTONIC, &tbegin);

  for (t = 0; t < num_stencil_iters; t++) {
    for (i = 1; i < n - 1; i++) {
      JBI1D_STENCIL(l2, l1);
    }
    l2[0] = l1[0];
    l2[n - 1] = l1[n - 1];
    memcpy(l1, l2, n * sizeof(*l2));
  }

  clock_gettime( CLOCK_MONOTONIC, &tend);
  for (int i = 0; i < n; i++) {
    jbi_out[i] = l1[i];
  }

  free(l1);
  free(l2);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
}

void jbi_init(double *jbi_in, double * jbi_out, int n) {
  int j;

  for (j = 0; j < n; j++) {
    jbi_in[j] = fabs(cos((double) j )) * (1 << 8);
    jbi_out[j] = 0;
  }
  jbi_in[0] = 0.0;
  jbi_in[n-1] = 0.0;
}


/* --------*/
/*  Tests  */
/* --------*/

void
djbi1d_half_diamonds_test(struct args_dimt args, double * jbi_in,
  double * jbi_out, struct benchscore * bsc)
{
  int pb_size = args.width, num_stencil_iters = args.iters;
  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  djbi1d_half_diamonds(pb_size, num_stencil_iters, jbi_in);
  clock_gettime( CLOCK_MONOTONIC, &tend);
  bsc->wallclock = ELAPSED_TIME(tend, tbegin);
  memcpy(jbi_out, jbi_in, sizeof(double) * pb_size);
}


void
djbi1d_sk_full_tiles_test(struct args_dimt args, double * jbi_in,
  double * jbi_out, struct benchscore * bsc)
{
  int pb_size = args.width, num_stencil_iters = args.iters;
  int i;
  int num_steps = (num_stencil_iters / T_ITERS) + 1;
  int num_strips = (pb_size/(T_WIDTH_DBL)) + 1;

  #ifdef DEBUG
    printf("iters : %i, n : %i -- %i num_steps, %i num_strips\n", num_steps,
      num_strips, num_stencil_iters, pb_size);
  #endif

  double * jbi_dashs = (double*) malloc(sizeof(double) *
    ((num_strips + num_steps) * T_WIDTH_DBL));
  double * jbi_slashs = (double*) malloc(sizeof(double) * num_steps * 2 * T_ITERS);



  /* Unused except in omp task pragmas. Attribute set to remove compiler
   * warnings
   */
  uint8_t ** tasks __attribute__ ((unused));

  if (jbi_dashs == NULL || jbi_slashs == NULL) {
    fprintf(stderr, "Error while allocating 2D arrays for skewed_tiles\n");
    return;
  }

  for (i = 0; i < pb_size; i ++) {
    jbi_dashs[i] = jbi_in[i];
  }

  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  tasks = djbi1d_sk_full_tiles(num_strips, num_steps, jbi_dashs, jbi_slashs);
  clock_gettime( CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);

  #ifdef DEBUG
    int t;
    if (task_index(tasks, num_strips, num_steps) > 0) {
      printf("The task index for dependencies doesn't seem correct ...\n");
      for (t = 0; t < num_steps; t ++) {
        for (i = 1; i < num_strips; i ++) {
          printf("%2i", tasks[t][i]);
        }
        printf("\n");
      }
    }
  #endif
  int start_stripe_top = num_steps * T_WIDTH_DBL;
  for (i = 0 ; i < pb_size; i ++) {
    jbi_out[i] = jbi_dashs[i + start_stripe_top];
  }

  free(jbi_dashs);
  free(jbi_slashs);
}

void
djbi1d_skewed_tiles_test(struct args_dimt args, double * jbi_in,
  double * jbi_out, struct benchscore * bsc)
{
  int pb_size = args.width, num_stencil_iters = args.iters;
  int i;
  int num_steps = (num_stencil_iters / T_ITERS) + 1;
  int num_strips = (pb_size/(T_WIDTH_DBL)) + 1;

  #ifdef DEBUG
    printf("iters : %i, n : %i -- %i num_steps, %i num_strips\n",
      num_steps, num_strips, num_stencil_iters, pb_size);
  #endif

  double * jbi_dashs = (double*) malloc(sizeof(double) *
                          ((num_strips + num_steps) * T_WIDTH_DBL));
  double * jbi_slashs = (double*) malloc(sizeof(double) *
                          num_steps * 2 * T_ITERS);

  if (jbi_dashs == NULL || jbi_slashs == NULL) {
    fprintf(stderr, "Error while allocating 2D arrays for skewed_tiles\n");
    return;
  }

  for (i = 0; i < pb_size; i ++) {
    jbi_dashs[i] = jbi_in[i];
  }

  clock_gettime( CLOCK_MONOTONIC, &tbegin);
  djbi1d_skewed_tiles(num_strips, num_steps, jbi_dashs, jbi_slashs);
  clock_gettime( CLOCK_MONOTONIC, &tend);

  bsc->wallclock = ELAPSED_TIME(tend, tbegin);

  int start_stripe_top = num_steps * T_WIDTH_DBL;
  for (i = 0 ; i < pb_size; i ++) {
    jbi_out[i] = jbi_dashs[i + start_stripe_top];
  }

  free(jbi_dashs);
  free(jbi_slashs);
}


/* Task functions :
 * do_i0_t0 : starting task, botile_tom-left corner of the space-time matrix
 * do_i0_t : left column, triangular tile
 * do_i_t0 : first time iteration, botile_tom line
 * do_i_t : regular task
 * do_in_t : right column, triangular tile
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

  for (i = 1; i < T_WIDTH_DBL + 1; i++) {
    l1[i] = dashs[i-1];
  }

  for (t = 0; t < T_ITERS * 2; t+=2) {
    l1[0] = 0.0;
    int right = max(T_WIDTH_DBL - t/2, 1);

    memcpy(l2, l1,(T_WIDTH_DBL + 1)  * sizeof(*l1));
    for (i = 1; i < right; i++) {
      l2[i] = (l1[i -1] + l1[i] + l1[i+1]) / 3.0 ;
    }

    slashs[t] = l1[right - 1];
    slashs[t + 1] = l1[right];
    memcpy(l1, l2,(T_WIDTH_DBL + 1)  * sizeof(*l1));
  }
  #ifdef DEBUG
    for (i = 0; i < DISPLAY_SIZE; i++)
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

  for (i = 1; i < T_WIDTH_DBL + 1; i++) {
   l1[i] = dashs[strpno * T_WIDTH_DBL + i - 1];
  }

  for (t = 0; t < 2 * T_ITERS; t += 2) {
    l1[0] = 0.0;
    int right = max(T_WIDTH_DBL - t/2, 0);

    memcpy(l2, l1,(T_WIDTH_DBL + 1)  * sizeof(*l1));
    for (i = 1; i < right; i++) {
      l2[i] = ((l1[i - 1] + l1[i] + l1[i + 1]) / 3.0);
    }
    slashs[tile_t * 2 * T_ITERS + t] = l1[right - 1];
    slashs[tile_t * 2 * T_ITERS + t + 1] = l1[right];

    memcpy(l1, l2,(T_WIDTH_DBL + 1)  * sizeof(*l1));
  }


  free(l1);
  free(l2);
}

inline void do_i_t(double * dashs, double * slashs, int strpno, int tile_t) {

  double * l1 = alloc_line(T_WIDTH_DBL + 2);
  double * l2 = alloc_line(T_WIDTH_DBL + 2);

  uint8_t t,i;

  #ifdef DEBUG
    if (tile_t == ((DBG_ITER / T_ITERS ))) {
      printf("Final line, task %i %i\n", tile_t, strpno - tile_t);
    }
  #endif

  /* Load dash */
  for (i = 2; i < T_WIDTH_DBL + 2; i++) {
    l1[i] = dashs[strpno * T_WIDTH_DBL + i - 2];
  }

  for (t = 0; t < T_ITERS * 2; t+=2) {
    /* Load slash part */
    l1[0] = slashs[tile_t * 2 * T_ITERS + t];
    l1[1] = slashs[tile_t * 2 * T_ITERS + t + 1];

    memcpy(l2, l1,(T_WIDTH_DBL + 2)  * sizeof(*l1));
    for (i = 2; i < T_WIDTH_DBL + 2; i++) {
      l2[i] = ((l1[i -2] + l1[i - 1] + l1[i]) / 3.0) ;
    }

    /* Write slash */
    slashs[tile_t * 2 * T_ITERS + t] = l1[T_WIDTH_DBL];
    slashs[tile_t * 2 * T_ITERS + t + 1] = l1[T_WIDTH_DBL + 1];

    memcpy(l1, l2,(T_WIDTH_DBL + 2)  * sizeof(*l1));
  }
  /* Write dash */
  for (i = 2; i < T_WIDTH_DBL; i++) {
    dashs[strpno * T_WIDTH_DBL + i - 2] = l2[i];
  }


  free(l1);
  free(l2);
}


inline void do_in_t(double * dashs, double * slashs, int strpno, int tile_t) {

  double * l1 = alloc_line(T_WIDTH_DBL + 2);
  double * l2 = alloc_line(T_WIDTH_DBL + 2);

  uint8_t t,i;

  #ifdef DEBUG
    if (tile_t == ((DBG_SIZE / T_WIDTH_DBL ))) {
      printf("Final task %i %i\n", tile_t, strpno - tile_t);
    } else {
      printf("Do do_in_t %i %i\n", tile_t, strpno - tile_t);
    }
  #endif

  for (t = 0; t < T_ITERS * 2; t+=2) {
    /* Load slash part */
    l1[0] = slashs[tile_t * 2 * T_ITERS + t];
    l1[1] = slashs[tile_t * 2 * T_ITERS + t + 1];
    int r = (2 + t/2);
    l1[r] = 0.0;
    for (i = 2; i <= r; i++) {
      l2[i] = ((l1[i -2] + l1[i - 1] + l1[i]) / 3.0) ;
    }
    memcpy(l1, l2,(T_WIDTH_DBL + 2)  * sizeof(*l1));
  }
  /* Write dash */
  for (i = 2; i < T_WIDTH_DBL; i++) {
    dashs[strpno * T_WIDTH_DBL + i - 2] = l1[i];
  }

  free(l1);
  free(l2);
}
