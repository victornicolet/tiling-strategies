#define _POSIX_C_SOURCE 200112L

#include <immintrin.h>
#include <inttypes.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "jacobi1d.h"

static struct timespec tend;
static struct timespec tbegin;

/* Functions describing different tasks in the computation
  Always inlined in the main body */
static inline void do_i_t(double *, double *, int, int)
    __attribute__((always_inline));
static inline void do_i0_t(double *, double *, int, int)
    __attribute__((always_inline));
static inline void do_i0_t0(double *, double *, int, int)
    __attribute__((always_inline));
static inline void do_in_t(double *, double *, int, int)
    __attribute__((always_inline));
static inline void do_top_hdiam(int, int, int, double **, double *)
    __attribute__((always_inline));
static inline void do_base_hdiam(int, int, int, double *, double **)
    __attribute__((always_inline));
static inline void do_topleft_hdiam(int, double **, double *)
    __attribute__((always_inline));
static inline void do_topright_hdiam(int, int, double **, double *)
    __attribute__((always_inline));

/*
* Different versions of jacobi1d have been implemented here.
*   - with half - diamonds : useful if we have few iterations, it allows
*      parallelism and a concurrent start.
*   - skewed tiles, with tasks
*   - with overlapping tiles
*   - sequential with swapping
*/

void djbi1d_hdiam_tasked(int pb_size, int num_iters, double *jbi,
                         double *jbi_out) {
  /* Tile bounds */
  int tile_no;
  int tile_base_sz = 2 * num_iters;
  int num_tiles = (pb_size / tile_base_sz);
  /* Store the border between base-down pyramids and base-up pyramids */
  double **tmp = alloc_double_mx(2, pb_size * sizeof(*tmp));

  uint8_t task_index[num_tiles] __attribute__((unused));

#pragma omp parallel
#pragma omp single
  {
/*
* Execute top-left task upfront : here the use of task pragma is only meant
* to satisfy the dependency.
*/

#pragma omp task depend(out : task_index[0])
    { do_base_hdiam(tile_no, num_iters, pb_size, jbi, tmp); }
    do_topleft_hdiam(num_iters, tmp, jbi_out);

    /* Loop over all tasks */
    for (tile_no = 1; tile_no < num_tiles; tile_no++) {
/* Base down tile */
#pragma omp task firstprivate(tile_no) shared(tmp) \
    depend(out : task_index[tile_no])
      { do_base_hdiam(tile_no, num_iters, pb_size, jbi, tmp); }

/* Tip down tile */
#pragma omp task firstprivate(tile_no) shared(tmp) \
    depend(in : task_index[tile_no - 1], task_index[tile_no])
      { do_top_hdiam(tile_no, num_iters, pb_size, tmp, jbi_out); }
    }

    /* Final task (top-right)*/

    do_top_hdiam(num_tiles, num_iters, pb_size, tmp, jbi_out);

  }
  /* End parallel - single region */
}

void djbi1d_hdiam_grouped(int pb_size, int num_iters, int num_procs,
                          double *jbi, double *jbi_out) {
  int tile_no, grp_no;
  num_procs = 2 * num_procs;
  /* Tile bounds */
  int i, t;
  int x0, r0, l0, r, l;
  double *li1, *li0;
  int tile_max;
  int tile_base_sz = 2 * num_iters;
  int num_tiles = (pb_size / tile_base_sz);
  /* Gourp size */
  int group_size = GROUP_FACTOR * (num_procs * L1_CACHE_SIZE) /
    (num_iters * sizeof(double) * 4) ;
  int num_grps = (num_tiles - 1) / group_size;
  /* Store the border between base-down pyramids and base-up pyramids */
  double **tmp = alloc_double_mx(2, pb_size * sizeof(*tmp));

  /* First execute first base-down tile and top-left corner */
  /* First base-down tile */
  do_base_hdiam(0, num_iters, pb_size, jbi, tmp);
  do_topleft_hdiam(num_iters, tmp, jbi_out);

  for (grp_no = 0; grp_no < num_grps + 1; grp_no++) {
    tile_max = min((grp_no + 1) * group_size + 1, num_tiles);

#pragma omp parallel for schedule(static) shared(tmp) \
    private(tile_no, x0, r0, l0, r, l, i ,t, li0, li1) \
    firstprivate(grp_no)

    for (tile_no = grp_no * group_size + 1; tile_no < tile_max; tile_no++) {
      /* Initial values */
      l0 = max(tile_no * tile_base_sz, 0);
      r0 = min(l0 + tile_base_sz, pb_size - 1);

      li1 = tmp[0] + l0;
      li0 = tmp[1] + l0;

      for (i = l0; i < r0; i++) {
        li0[i - l0] = jbi[i];
        li1[i - l0] = 0.0;
      }

      for (t = 0; t < num_iters - 1; t++) {
        l = max(l0 + t + 1, 1);
        r = min(l0 + tile_base_sz - t - 1, pb_size);

        for (i = l; i < r; i++) {
          li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
        }

        swap(&li0, &li1);
      }
    }

#pragma omp parallel for schedule(static) shared(tmp) \
    private(tile_no, x0, r0, l0, r, l, i ,t, li0, li1) \
    firstprivate(grp_no)

    for (tile_no = grp_no * group_size + 1; tile_no < tile_max; tile_no++) {
      x0 = tile_no * tile_base_sz;
      l0 = max(x0 - num_iters - 1, 0);
      r0 = min(x0 + num_iters + 1, pb_size - 1);

      li1 = tmp[0] + l0;
      li0 = tmp[1] + l0;

      for (t = 0; t < num_iters; t++) {
        l = max(x0 - (t + 1), 1);
        r = min(x0 + t, pb_size - 2);
        /* Load from the border-storing array */

        for (i = l; i <= r; i++) {
          li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
        }

        swap(&li0, &li1);
      }

      /* Copy back to memory */
      for (i = l0 + 1; i < r0 - 1; i++) {
        jbi_out[i] = li0[i - l0];
      }
    }
  }


  free_mx((void **)tmp, 2);
}

void djbi1d_half_diamonds(int pb_size, int num_iters, double *jbi_in,
                          double *jbi_out) {
  int i, t;
  int x0, r0, l0, r, l;
  int tile_base_sz = num_iters * 2;
  double *li1, *li0;
  int tile_no;
  int num_tiles = (pb_size / (2 * num_iters));
  /* Store the border between base-down pyramids and base-up pyramids */
  double **tmp = alloc_double_mx(2, pb_size * sizeof(*tmp));

  /* First loop : base down tiles */
#pragma omp parallel for schedule(static) shared(tmp) \
  private(tile_no, x0, r0, l0, r, l, i ,t, li0, li1)
  for (tile_no = 0; tile_no < num_tiles; tile_no++) {
    /* Initial values */
    l0 = max(tile_no * tile_base_sz, 0);
    r0 = min(l0 + tile_base_sz, pb_size - 1);

    li1 = tmp[0] + l0;
    li0 = tmp[1] + l0;

    for (i = l0; i < r0; i++) {
      li0[i - l0] = jbi_in[i];
      li1[i - l0] = 0.0;
    }

    for (t = 0; t < num_iters - 1; t++) {
      l = max(l0 + t + 1, 1);
      r = min(l0 + tile_base_sz - t - 1, pb_size);

      for (i = l; i < r; i++) {
        li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
      }

      swap(&li0, &li1);
    }
  }

  /* Half-tile at beggining, in front of parallel loop */
 do_topleft_hdiam(num_iters, tmp, jbi_out);


  /* Second loop : tip down tiles */
#pragma omp parallel for schedule(static) shared(tmp) \
  private(tile_no, x0, r0, l0, r, l, i ,t, li0, li1)
  for (tile_no = 1; tile_no < num_tiles; tile_no++) {
      x0 = tile_no * tile_base_sz;
      l0 = max(x0 - num_iters - 1, 0);
      r0 = min(x0 + num_iters + 1, pb_size - 1);

      li1 = tmp[0] + l0;
      li0 = tmp[1] + l0;

      for (t = 0; t < num_iters; t++) {
        l = max(x0 - (t + 1), 1);
        r = min(x0 + t, pb_size - 2);
        /* Load from the border-storing array */

        for (i = l; i <= r; i++) {
          li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
        }

        swap(&li0, &li1);
      }

      /* Copy back to memory */
      for (i = l0 + 1; i < r0 - 1; i++) {
        jbi_out[i] = li0[i - l0];
      }
  }

  do_topright_hdiam(num_iters, pb_size, tmp, jbi_out);
  free_mx((void **)tmp, 2);

}

void ljbi1d_half_diamonds(int pb_size, int num_iters, long *jbi,
                          long *jbi_out) {
  int tile_no, t, i;
  // Tile bounds
  int l, l0, r, r0, x0;

  int tile_base_sz = 2 * num_iters;
  int num_tiles = (pb_size / tile_base_sz);
  // Store the border between base-down pyramids and base-up pyramids
  long **tmp = alloc_long_mx(2, pb_size * sizeof(*tmp));

  for (i = 0; i < pb_size; i++) {
    tmp[1][i] = 9999;
    tmp[0][i] = 9999;
  }

/* First loop : base down tiles */

#pragma omp parallel for schedule(static) shared(tmp) private(l0, r0, l, r, x0)

  for (tile_no = 0; tile_no < num_tiles; tile_no++) {
    long li1[tile_base_sz], li0[tile_base_sz];

    /* Initial values */
    l0 = max(tile_no * tile_base_sz, 0);
    r0 = min(l0 + tile_base_sz, pb_size - 1);
    for (i = l0; i < r0; i++) {
      li0[i - l0] = jbi[i];
      li1[i - l0] = 0.0;
    }

    tmp[0][l0] = li0[0];
    tmp[1][l0 + 1] = li0[1];
    tmp[0][r0 - 1] = li0[r0 - l0 - 1];
    tmp[1][r0 - 2] = li0[r0 - l0 - 2];

    for (t = 0; t < num_iters - 1; t++) {
      l = max(l0 + t + 1, 1);
      r = min(l0 + tile_base_sz - t - 1, pb_size);

      for (i = l; i < r; i++) {
        li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
      }

      tmp[0][r - 1] = li1[r - l0 - 1];
      tmp[0][l] = li1[l - l0];
      if (r - 2 != l) {
        tmp[1][r - 2] = li1[r - l0 - 2];
        tmp[1][l + 1] = li1[l - l0 + 1];
      }

      memcpy(li0, li1, (tile_base_sz) * sizeof(*li1));
    }
  }

  /* Half-tile at beggining */
  long li1[tile_base_sz], li0[tile_base_sz];

  for (i = 0; i < tile_base_sz; i++) {
    li1[i] = 999;
  }

  li1[0] = tmp[0][0];
  li0[0] = tmp[0][0];
  li0[1] = tmp[1][1];
  for (t = 0; t < num_iters; t++) {
    /* Load from the border-storing array */
    li0[t + 1] = tmp[1][t + 1];
    li0[t] = tmp[0][t];

    for (i = 1; i < t + 1; i++) {
      li1[i] = (li0[i - 1] + li0[i] + li0[i + 1]) / 3.0;
    }

    memcpy(li0, li1, (tile_base_sz) * sizeof(*li1));
  }
  for (i = 0; i < num_iters + 1; i++) {
    jbi_out[i] = li0[i];
  }

/* Second loop : tip down tiles */

#pragma omp parallel for schedule(static) shared(tmp) private(l0, r0, l, r, x0)

  for (tile_no = 1; tile_no < num_tiles + 1; tile_no++) {
    x0 = tile_no * tile_base_sz;
    l0 = max(x0 - num_iters - 1, 0);
    r0 = min(x0 + num_iters + 1, pb_size);

    long li1[tile_base_sz + 2], li0[tile_base_sz + 2];

    for (t = 0; t < num_iters; t++) {
      l = max(x0 - (t + 1), 1);
      r = min(x0 + t, pb_size - 1);
      /* Load from the border-storing array */
      li0[l - l0 - 1] = tmp[1][l - 1];
      li0[r - l0 + 1] = tmp[1][r + 1];
      li0[r - l0] = tmp[0][r];
      li0[l - l0] = tmp[0][l];

      for (i = l; i <= r; i++) {
        li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
      }
      memcpy(li0, li1, (tile_base_sz + 2) * sizeof(*li1));
    }
    /* Copy back to memory */
    for (i = l0 + 1; i < r0; i++) {
      jbi_out[i] = li0[i - l0];
    }
  }
  free(tmp);
}

/* In this algorithm we work only on full parallelograms : no partial tiles.
 * We are interested in performance measures, and understanding, rather
 * than corectness here
 */
uint8_t **djbi1d_sk_full_tiles(int num_strips, int num_steps, double *dashs,
                               double *slashs) {
  int i, t;

  uint8_t **taskdep_index = malloc(num_steps * sizeof(*taskdep_index));
  for (i = 0; i < num_strips; i++) {
    taskdep_index[i] = malloc(num_steps * sizeof *taskdep_index[i]);
  }

#ifndef SEQ
#pragma omp parallel
#endif
  {
#ifndef SEQ
#pragma omp single
#endif
    {
      for (t = 0; t < num_steps; t++) {
        for (i = 1; i < num_strips; i++) {
          taskdep_index[t][i] = 0;
          int strpno = i + t;
          if (t == 0 && i == 1) {
#ifndef SEQ
#pragma omp task firstprivate(i, t) depend(out : taskdep_index[t][i])
#endif
            {
              do_i_t(dashs, slashs, 1, 0);
              taskdep_index[t][i] ^= 1;
            }
          } else if (t == 0 && i > 1) {
#ifndef SEQ
#pragma omp task firstprivate(i, t) depend(out : taskdep_index[t][i]) \
                                  depend(in : taskdep_index[i - 1][0])
#endif
            {
              do_i_t(dashs, slashs, strpno, 0);
              taskdep_index[t][i] ^= 1;
            }
          } else {
#ifndef SEQ
#pragma omp task firstprivate(i, t) depend(out : taskdep_index[t][i]) depend( \
    in : taskdep_index[t][i - 1], taskdep_index[t - 1][i + 1])
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

void djbi1d_skewed_tiles(int num_strips, int num_steps, double *dashs,
                         double *slashs) {
  /* In this version we assume T_ITERS = T_WIDTH_DBL so we have to make a
   *   difference between regular tiles ( parallelogram-shaped ones) and
   *   triangular tiles, but the patile_tern is quite regular and dependencies
   *are
   *   straightforward at the boundaries of the domain.
   *   ---------------
   *
   *  num_steps = (iters / T_ITERS) + 1;
   *   num_strips = (n / T_WIDTH_DBL) + 1 ;
   */

  /* Unused except in omp task pragmas */
  uint8_t taskdep_index[num_strips + 1][num_steps] __attribute__((unused));

#ifdef DEBUG_PARALLEL
  int *tsk = (int *)malloc(sizeof(int) * (num_strips + 1) * num_steps);
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

        if (tile_t == 0 && tile_i == 0) {
#ifndef SEQ
#pragma omp task firstprivate(tile_i, tile_t) depend(out : taskdep_index[0][0])
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
#pragma omp task firstprivate(tile_i, tile_t) depend( \
    in : taskdep_index[tile_i - 1][0])                \
        depend(out : taskdep_index[tile_i][tile_t])
#endif
          {
#ifdef DEBUG_PARALLEL
            if (tsk[tile_i - 1] != 1) {
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
#pragma omp task firstprivate(tile_i, tile_t) depend( \
    in : taskdep_index[1][tile_t - 1])                \
        depend(out : taskdep_index[tile_i][tile_t])
#endif
          {
#ifdef DEBUG_PARALLEL
            if (tsk[(tile_t - 1) * num_steps + tile_i] != 1) {
              printf("Unsatisified dependency !\n");
            }
            tsk[tile_t * num_steps + tile_i] = 1;
#endif
            do_i0_t(dashs, slashs, strpno, tile_t);
          }
        } else if (tile_i == num_strips) {
#ifndef SEQ
#pragma omp task firstprivate(tile_i, tile_t) depend( \
    in : taskdep_index[tile_i - 1][tile_t])           \
        depend(out : taskdep_index[tile_i][tile_t])
#endif
          {
#ifdef DEBUG_PARALLEL
            if (tsk[tile_t * num_steps + tile_i - 1] != 1) {
              printf("Unsatisified dependency !\n");
            }
            tsk[tile_t * num_steps + tile_i] = 1;
#endif
            do_in_t(dashs, slashs, strpno, tile_t);
          }

        } else {
/* Regular tile two in and out dependencies */
#ifndef SEQ
#pragma omp task firstprivate(tile_i, tile_t) depend(                    \
    in : taskdep_index[tile_i - 1][tile_t],                              \
                                   taskdep_index[tile_i + 1][tile_t -    \
                                                             1]) depend( \
                                       out : taskdep_index[tile_i][tile_t])
#endif
          {
#ifdef DEBUG_PARALLEL
            if (tsk[tile_t * num_steps + tile_i - 1] != 1 ||
                tsk[(tile_t - 1) * num_steps + (tile_i + 1)] != 1) {
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

double djbi1d_diamond_tiles(int n, int num_stencil_iters, double **jbi) {
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
        l = max(tile_i * T_WIDTH_DBL_DIAM - (t - bot), 1);
        r = min((tile_i + 1) * T_WIDTH_DBL_DIAM - (t - bot), n - 1);
        for (int i = l; i < r; i++) {
          JBI1D_STENCIL_T(jbi);
        }
      }
    }
  }

  /* Not implemented */
  return -1.0;
}

double djbi1d_omp_naive(struct args_dimt args, double *jbi_in,
	double *jbi_out) {
  int n = args.width;
  int num_stencil_iters = args.iters;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  int t, i;
  double *l1 = (double *)aligned_alloc(CACHE_LINE_SIZE, sizeof(double) * n);
  double *l2 = (double *)aligned_alloc(CACHE_LINE_SIZE, sizeof(double) * n);
  memcpy(l1, jbi_in, n * sizeof(double));

  clock_gettime(CLOCK_MONOTONIC, &tbegin);

  for (t = 0; t < num_stencil_iters; t++) {
#ifdef SEQ
#pragma omp parallel for schedule(static)
#endif
    for (i = 1; i < n - 1; i++) {
      JBI1D_STENCIL(l2, l1);
    }
    swap(&l1, &l2);
  }

  clock_gettime(CLOCK_MONOTONIC, &tend);

  for (int i = 0; i < n; i++) {
    jbi_out[i] = l1[i];
  }

  free(l1);
  free(l2);

  elapsed = ELAPSED_TIME(tend, tbegin);

  return elapsed;
}

double djbi1d_omp_overlap(struct args_dimt args, double *jbi_in,
  double *jbi_out) {
  int pb_size = args.width, num_stencil_iters = args.iters;
  int tile_i, tile_t, t, i;
  int tile_base_sz = T_WIDTH_DBL_OVERLAP + T_ITERS * 2;

  clock_gettime(CLOCK_MONOTONIC, &tbegin);

  for (tile_t = 0; tile_t <= num_stencil_iters / T_ITERS; tile_t++) {
#ifndef SEQ
#pragma omp parallel for schedule(static)
#endif
    for (tile_i = 0; tile_i <= (pb_size / T_WIDTH_DBL_OVERLAP); tile_i++) {
      double *lvl1 = malloc(tile_base_sz * sizeof(*lvl1));
      double *lvl0 = malloc(tile_base_sz * sizeof(*lvl0));

      /* Compute tile bounds */
      int bot = max((T_ITERS * tile_t), 1);
      int top = min((T_ITERS * (tile_t + 1)), num_stencil_iters);
      int h = top - bot;
      int l0 = max((T_WIDTH_DBL_OVERLAP * tile_i), 0);
      int r0 = min((T_WIDTH_DBL_OVERLAP * (tile_i + 1)), pb_size);
      int l = max((l0 - h), 0);
      int r = min((r0 + h), pb_size);
      {
        /* Read tile base */
        for (i = l; i < r; i++) {
          lvl0[i - l] = jbi_in[i];
          lvl1[i - l] = 0.0;
        }

        for (t = 0; t < h; t++) {
          int lt = max(l0 - t, 1);
          int rt = min(r0 + t, pb_size);

          for (i = lt; i < rt; i++) {
            lvl1[i - l] =
                (lvl0[i - l - 1] + lvl0[i - l] + lvl0[i - l + 1]) / 3.0;
          }
          swap(&lvl0, &lvl1);
        }

        /* Write tile top */
        for (i = l0; i < r0; i++) {
          jbi_out[i] = lvl0[i - l];
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

  clock_gettime(CLOCK_MONOTONIC, &tend);

  return ELAPSED_TIME(tend, tbegin);
}

double djbi1d_sequential(struct args_dimt args, double *jbi_in,
  double *jbi_out) {
  int num_stencil_iters = args.iters, n = args.width;
  /* Boundaries initial condition */
  int t, i;
  double *l1 = (double *)aligned_alloc(CACHE_LINE_SIZE, sizeof(*l1) * n);
  double *l2 = (double *)aligned_alloc(CACHE_LINE_SIZE, sizeof(*l2) * n);
  memcpy(l1, jbi_in, n * sizeof(*jbi_in));

  clock_gettime(CLOCK_MONOTONIC, &tbegin);

  for (t = 0; t < num_stencil_iters; t++) {
    for (i = 1; i < n - 1; i++) {
      JBI1D_STENCIL(l2, l1);
    }
    l2[0] = l1[0];
    l2[n - 1] = l1[n - 1];
    swap(&l1, &l2);
  }

  clock_gettime(CLOCK_MONOTONIC, &tend);
  for (int i = 0; i < n; i++) {
    jbi_out[i] = l1[i];
  }

  free(l1);
  free(l2);

  return ELAPSED_TIME(tend, tbegin);

}

double ljbi1d_sequential(struct args_dimt args, long *jbi_in, long *jbi_out) {
  int num_stencil_iters = args.iters, n = args.width;
  /* Boundaries initial condition */
  int t, i;
  long *l1, *l2;

  l1 = aligned_alloc(CACHE_LINE_SIZE, sizeof(*l1) * n);
  l2 = aligned_alloc(CACHE_LINE_SIZE, sizeof(*l2) * n);
  memcpy(l1, jbi_in, n * sizeof(*jbi_in));

  clock_gettime(CLOCK_MONOTONIC, &tbegin);

  for (t = 0; t < num_stencil_iters; t++) {
    for (i = 1; i < n - 1; i++) {
      JBI1D_STENCIL(l2, l1);
    }
    l2[0] = l1[0];
    l2[n - 1] = l1[n - 1];
    memcpy(l1, l2, n * sizeof(*l2));
  }

  clock_gettime(CLOCK_MONOTONIC, &tend);
  for (int i = 0; i < n; i++) {
    jbi_out[i] = l1[i];
  }

  free(l1);
  free(l2);

  return ELAPSED_TIME(tend, tbegin);
}

double djbi1d_from_pluto(struct args_dimt args, double *jbi_in,
  double *jbi_out) {

  int N,T;
  int c1, c2, c3, c4, c5;
  int i, j, k, l, t;
  (void)i;
  (void)j;
  (void)k;
  (void)l;
  (void)t;
  register int lb, ub;

  N = args.width;
  T = args.iters;

  double *a, *b;
  a = jbi_in;
  b = jbi_out;

  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  /* Generated from jacobi-imper.sched.cloog by CLooG v0.14.1 64 bits in 0.01s. */
 for (c1=-1;c1<=floord(N+3*T-4,2048);c1++) {
     lb = max(max(ceild(2048*c1-T+1,2048),ceild(4096*c1-2045,6144)),0);
     ub = min(min(floord(2048*c1+2047,2048),floord(4096*c1+N+4093,6144)),floord(N+2*T-3,2048));

#pragma omp parallel for shared(c1,lb,ub,a,b) private(c2,c3,c4,c5,i,j,k,l)
   for (c2=lb;c2<=ub;c2++) {
     if (c1 >= max(c2,ceild(6144*c2-N+2,4096))) {
       c3 = 4096*c1-4096*c2 ;
        for (c4=max(4096*c1-4096*c2+2,2048*c2);c4<=min(4096*c1-4096*c2+N-2,2048*c2+2047);c4++) {
          c5 = 0 ;
          if ((c1-c2)%2 == 0) {
            i = (c1-c2)/2 ;
            j = -c1+2*c2 ;
            k = 2048*c1-2048*c2 ;
            l = -4096*c1+4096*c2+c4 ;
            S1((c1-c2)/2,-c1+2*c2,2048*c1-2048*c2,-4096*c1+4096*c2+c4) ;
          }
        }
      }
      if ((c1 <= floord(4096*c2-1,4096)) && (c2 <= floord(N-2,2048))) {
        c3 = 0 ;
        for (c4=max(2048*c2,2);c4<=min(2048*c2+2047,N-2);c4++) {
          c5 = 0 ;
          if ((c1-c2)%2 == 0) {
            i = (c1-c2)/2 ;
            j = -c1+2*c2 ;
            k = 0 ;
            l = c4 ;
            S1((c1-c2)/2,-c1+2*c2,0,c4) ;
          }
        }
      }
      for (c3=max(max(1,2048*c2-N+2),4096*c1-4096*c2+1);c3<=min(min(4096*c1-4096*c2+4094,2048*c2+2045),2*T-2);c3++) {
        for (c4=max(2048*c2,c3+2);c4<=min(c3+N-2,2048*c2+2047);c4++) {
          c5 = 0 ;
          if ((c1-c2)%2 == 0) {
            i = (c1-c2)/2 ;
            j = -c1+2*c2 ;
            if (c3%2 == 0) {
              k = c3/2 ;
              l = -c3+c4 ;
              S1((c1-c2)/2,-c1+2*c2,c3/2,-c3+c4) ;
            }
          }
          c5 = 1 ;
          if ((c1-c2)%2 == 0) {
            i = (c1-c2)/2 ;
            j = -c1+2*c2 ;
            if ((c3-1)%2 == 0) {
              k = (c3-1)/2 ;
              l = -c3+c4 ;
              S2((c1-c2)/2,-c1+2*c2,(c3-1)/2,-c3+c4) ;
            }
          }
        }
      }
      if (c1 <= min(floord(3072*c2-1025,2048),floord(2048*c2+T-2048,2048))) {
        c3 = 4096*c1-4096*c2+4095 ;
        for (c4=max(4096*c1-4096*c2+4097,2048*c2);c4<=min(4096*c1-4096*c2+N+4093,2048*c2+2047);c4++) {
          c5 = 1 ;
          if ((c1-c2)%2 == 0) {
            i = (c1-c2)/2 ;
            j = -c1+2*c2 ;
            k = 2048*c1-2048*c2+2047 ;
            l = -4096*c1+4096*c2+c4-4095 ;
            S2((c1-c2)/2,-c1+2*c2,2048*c1-2048*c2+2047,-4096*c1+4096*c2+c4-4095) ;
          }
        }
      }
      if ((c1 >= ceild(4096*c2+2*T-4095,4096)) && (c2 >= ceild(T-1023,1024))) {
        c3 = 2*T-1 ;
        for (c4=max(2*T+1,2048*c2);c4<=min(2048*c2+2047,N+2*T-3);c4++) {
          c5 = 1 ;
          if ((c1-c2)%2 == 0) {
            i = (c1-c2)/2 ;
            j = -c1+2*c2 ;
            k = T-1 ;
            l = c4-2*T+1 ;
            S2((c1-c2)/2,-c1+2*c2,T-1,c4-2*T+1) ;
          }
        }
      }
    }
  }
  clock_gettime(CLOCK_MONOTONIC, &tend);

  return ELAPSED_TIME(tend, tbegin);
}

/* ==========================================================================*/
/*                               Tests                                       */
/* ==========================================================================*/

double djbi1d_hdiam_tasked_test(struct args_dimt args, double *jbi_in,
  double *jbi_out) {
  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  djbi1d_hdiam_tasked(args.width, args.iters, jbi_in, jbi_out);
  clock_gettime(CLOCK_MONOTONIC, &tend);

  return ELAPSED_TIME(tend, tbegin);
}

double djbi1d_half_diamonds_test(struct args_dimt args, double *jbi_in,
  double *jbi_out) {
  int pb_size = args.width, num_stencil_iters = args.iters;
  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  djbi1d_half_diamonds(pb_size, num_stencil_iters, jbi_in, jbi_out);
  clock_gettime(CLOCK_MONOTONIC, &tend);
  return ELAPSED_TIME(tend, tbegin);
}

double djbi1d_hdiam_grouped_test(struct args_dimt args, double *jbi_in,
  double *jbi_out) {
  int pb_size = args.width, num_stencil_iters = args.iters;
  int num_procs = DEFAULT_PROC_NUM;
#ifdef _SC_NPROCESSORS_ONLN
  num_procs = sysconf(_SC_NPROCESSORS_ONLN);
#endif
  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  djbi1d_hdiam_grouped(pb_size, num_stencil_iters, num_procs, jbi_in, jbi_out);
  clock_gettime(CLOCK_MONOTONIC, &tend);
  return ELAPSED_TIME(tend, tbegin);
}

double ljbi1d_half_diamonds_test(struct args_dimt args, long *jbi_in,
  long *jbi_out) {
  int pb_size = args.width, num_stencil_iters = args.iters;
  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  ljbi1d_half_diamonds(pb_size, num_stencil_iters, jbi_in, jbi_out);
  clock_gettime(CLOCK_MONOTONIC, &tend);
  return ELAPSED_TIME(tend, tbegin);
}

double djbi1d_sk_full_tiles_test(struct args_dimt args, double *jbi_in,
                               double *jbi_out) {
  int pb_size = args.width, num_stencil_iters = args.iters;
  int i;
  int num_steps = (num_stencil_iters / T_ITERS) + 1;
  int num_strips = (pb_size / (T_WIDTH_DBL)) + 1;

#ifdef DEBUG
  printf("iters : %i, n : %i -- %i num_steps, %i num_strips\n", num_steps,
         num_strips, num_stencil_iters, pb_size);
#endif

  double *jbi_dashs = (double *)malloc(
      sizeof(double) * ((num_strips + num_steps) * T_WIDTH_DBL));
  double *jbi_slashs =
      (double *)malloc(sizeof(double) * num_steps * 2 * T_ITERS);

  /* Unused except in omp task pragmas. Attribute set to remove compiler
   * warnings
   */
  uint8_t **tasks __attribute__((unused));

  if (jbi_dashs == NULL || jbi_slashs == NULL) {
    fprintf(stderr, "Error while allocating 2D arrays for skewed_tiles\n");
    return -1.0;
  }

  for (i = 0; i < pb_size; i++) {
    jbi_dashs[i] = jbi_in[i];
  }

  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  tasks = djbi1d_sk_full_tiles(num_strips, num_steps, jbi_dashs, jbi_slashs);
  clock_gettime(CLOCK_MONOTONIC, &tend);


#ifdef DEBUG
  int t;
  if (task_index(tasks, num_strips, num_steps) > 0) {
    printf("The task index for dependencies doesn't seem correct ...\n");
    for (t = 0; t < num_steps; t++) {
      for (i = 1; i < num_strips; i++) {
        printf("%2i", tasks[t][i]);
      }
      printf("\n");
    }
  }
#endif
  int start_stripe_top = num_steps * T_WIDTH_DBL;
  for (i = 0; i < pb_size; i++) {
    jbi_out[i] = jbi_dashs[i + start_stripe_top];
  }

  free(jbi_dashs);
  free(jbi_slashs);


  return ELAPSED_TIME(tend, tbegin);
}

double djbi1d_skewed_tiles_test(struct args_dimt args, double *jbi_in,
  double *jbi_out) {
  int pb_size = args.width, num_stencil_iters = args.iters;
  int i;
  int num_steps = (num_stencil_iters / T_ITERS) + 1;
  int num_strips = (pb_size / (T_WIDTH_DBL)) + 1;

#ifdef DEBUG
  printf("iters : %i, n : %i -- %i num_steps, %i num_strips\n", num_steps,
         num_strips, num_stencil_iters, pb_size);
#endif

  double *jbi_dashs = (double *)malloc(
      sizeof(double) * ((num_strips + num_steps) * T_WIDTH_DBL));
  double *jbi_slashs =
      (double *)malloc(sizeof(double) * num_steps * 2 * T_ITERS);

  if (jbi_dashs == NULL || jbi_slashs == NULL) {
    fprintf(stderr, "Error while allocating 2D arrays for skewed_tiles\n");
    return -1.0;
  }

  for (i = 0; i < pb_size; i++) {
    jbi_dashs[i] = jbi_in[i];
  }

  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  djbi1d_skewed_tiles(num_strips, num_steps, jbi_dashs, jbi_slashs);
  clock_gettime(CLOCK_MONOTONIC, &tend);


  int start_stripe_top = num_steps * T_WIDTH_DBL;
  for (i = 0; i < pb_size; i++) {
    jbi_out[i] = jbi_dashs[i + start_stripe_top];
  }

  free(jbi_dashs);
  free(jbi_slashs);

  return ELAPSED_TIME(tend, tbegin);
}

/* Task functions :
 * do_i0_t0 : starting task, botile_tom-left corner of the space-time matrix
 * do_i0_t : left column, triangular tile
 * do_i_t0 : first time iteration, botile_tom line
 * do_i_t : regular task
 * do_in_t : right column, triangular tile
 For half-diamonds version :
  * do_base_hdiam : base-down triangular tiles
  * do_top_hdiam : tip-down triangular tiles
  * do_topleft_hdiam : top-left half-tile
 */

static inline void do_i0_t0(double *dashs, double *slashs, int tile_t,
                            int tile_i) {
  double *l1 = alloc_line(T_WIDTH_DBL + 1);
  double *l2 = alloc_line(T_WIDTH_DBL + 1);

  uint8_t t, i;

#ifdef DEBUG
  printf("First task T I, %i %i\n", tile_t, tile_i);
#endif

  for (i = 1; i < T_WIDTH_DBL + 1; i++) {
    l1[i] = dashs[i - 1];
  }

  for (t = 0; t < T_ITERS * 2; t += 2) {
    l1[0] = 0.0;
    int right = max(T_WIDTH_DBL - t / 2, 1);

    memcpy(l2, l1, (T_WIDTH_DBL + 1) * sizeof(*l1));
    for (i = 1; i < right; i++) {
      l2[i] = (l1[i - 1] + l1[i] + l1[i + 1]) / 3.0;
    }

    slashs[t] = l1[right - 1];
    slashs[t + 1] = l1[right];
   swap(&l1, &l2);
  }
#ifdef DEBUG
  for (i = 0; i < DISPLAY_SIZE; i++) printf("%10.3f", dashs[i]);
  printf("\n");
#endif

  free(l1);
  free(l2);
}

static inline void do_i0_t(double *dashs, double *slashs, int strpno,
                           int tile_t) {
  double *l1 = alloc_line(T_WIDTH_DBL + 1);
  double *l2 = alloc_line(T_WIDTH_DBL + 1);

  uint8_t t, i;

#ifdef DEBUG
  printf("Do do_i0_t %i %i\n", tile_t, strpno - tile_t);
#endif

  for (i = 1; i < T_WIDTH_DBL + 1; i++) {
    l1[i] = dashs[strpno * T_WIDTH_DBL + i - 1];
  }

  for (t = 0; t < 2 * T_ITERS; t += 2) {
    l1[0] = 0.0;
    int right = max(T_WIDTH_DBL - t / 2, 0);

    memcpy(l2, l1, (T_WIDTH_DBL + 1) * sizeof(*l1));
    for (i = 1; i < right; i++) {
      l2[i] = ((l1[i - 1] + l1[i] + l1[i + 1]) / 3.0);
    }
    slashs[tile_t * 2 * T_ITERS + t] = l1[right - 1];
    slashs[tile_t * 2 * T_ITERS + t + 1] = l1[right];

    swap(&l1, &l2);
  }

  free(l1);
  free(l2);
}

static inline void do_i_t(double *dashs, double *slashs, int strpno,
                          int tile_t) {
  double *l1 = alloc_line(T_WIDTH_DBL + 2);
  double *l2 = alloc_line(T_WIDTH_DBL + 2);

  uint8_t t, i;

#ifdef DEBUG
  if (tile_t == ((DEBUG_ITER / T_ITERS))) {
    printf("Final line, task %i %i\n", tile_t, strpno - tile_t);
  }
#endif

  /* Load dash */
  for (i = 2; i < T_WIDTH_DBL + 2; i++) {
    l1[i] = dashs[strpno * T_WIDTH_DBL + i - 2];
  }

  for (t = 0; t < T_ITERS * 2; t += 2) {
    /* Load slash part */
    l1[0] = slashs[tile_t * 2 * T_ITERS + t];
    l1[1] = slashs[tile_t * 2 * T_ITERS + t + 1];

    memcpy(l2, l1, (T_WIDTH_DBL + 2) * sizeof(*l1));
    for (i = 2; i < T_WIDTH_DBL + 2; i++) {
      l2[i] = ((l1[i - 2] + l1[i - 1] + l1[i]) / 3.0);
    }

    /* Write slash */
    slashs[tile_t * 2 * T_ITERS + t] = l1[T_WIDTH_DBL];
    slashs[tile_t * 2 * T_ITERS + t + 1] = l1[T_WIDTH_DBL + 1];

    swap(&l1, &l2);
  }
  /* Write dash */
  for (i = 2; i < T_WIDTH_DBL; i++) {
    dashs[strpno * T_WIDTH_DBL + i - 2] = l2[i];
  }

  free(l1);
  free(l2);
}

static inline void do_in_t(double *dashs, double *slashs, int strpno,
                           int tile_t) {
  double *l1 = alloc_line(T_WIDTH_DBL + 2);
  double *l2 = alloc_line(T_WIDTH_DBL + 2);

  uint8_t t, i;

#ifdef DEBUG
  if (tile_t == ((DEBUG_SIZE / T_WIDTH_DBL))) {
    printf("Final task %i %i\n", tile_t, strpno - tile_t);
  } else {
    printf("Do do_in_t %i %i\n", tile_t, strpno - tile_t);
  }
#endif

  for (t = 0; t < T_ITERS * 2; t += 2) {
    /* Load slash part */
    l1[0] = slashs[tile_t * 2 * T_ITERS + t];
    l1[1] = slashs[tile_t * 2 * T_ITERS + t + 1];
    int r = (2 + t / 2);
    l1[r] = 0.0;
    for (i = 2; i <= r; i++) {
      l2[i] = ((l1[i - 2] + l1[i - 1] + l1[i]) / 3.0);
    }
    swap(&l1, &l2);
  }
  /* Write dash */
  for (i = 2; i < T_WIDTH_DBL; i++) {
    dashs[strpno * T_WIDTH_DBL + i - 2] = l1[i];
  }

  free(l1);
  free(l2);
}

/* Tasks for half-diamond version */

static inline void do_base_hdiam(int tile_no, int num_iters, int pb_size,
                                 double *jbi, double **tmp) {
  int i, t;
  int r0, l0, r, l;
  int tile_base_sz = num_iters * 2;
  double *li1, *li0;


  /* Initial values */
  l0 = max(tile_no * tile_base_sz, 0);
  r0 = min(l0 + tile_base_sz, pb_size - 1);

  li1 = tmp[0] + l0;
  li0 = tmp[1] + l0;

  for (i = l0; i < r0; i++) {
    li0[i - l0] = jbi[i];
    li1[i - l0] = 0.0;
  }

  for (t = 0; t < num_iters - 1; t++) {
    l = max(l0 + t + 1, 1);
    r = min(l0 + tile_base_sz - t - 1, pb_size);

    for (i = l; i < r; i++) {
      li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
    }

    swap(&li0, &li1);
  }
}

static inline void do_top_hdiam(int tile_no, int num_iters, int pb_size,
                                double **tmp, double *jbi_out) {
  int t, i;
  int x0, l0, r0, l, r;
  int tile_base_sz = 2 * num_iters;

  x0 = tile_no * tile_base_sz;
  l0 = max(x0 - num_iters - 1, 0);
  r0 = min(x0 + num_iters + 1, pb_size - 1);

  double *li1, *li0;

  li1 = tmp[0] + l0;
  li0 = tmp[1] + l0;

  for (t = 0; t < num_iters; t++) {
    l = max(x0 - (t + 1), 1);
    r = min(x0 + t, pb_size - 2);
    /* Load from the border-storing array */

    for (i = l; i <= r; i++) {
      li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
    }

    swap(&li0, &li1);
  }

  /* Copy back to memory */
  for (i = l0 + 1; i < r0 - 1; i++) {
    jbi_out[i] = li0[i - l0];
  }
}

static inline void do_topleft_hdiam(int num_iters, double **tmp,
  double *jbi_out) {
  int i, t;
  double *li1, *li0;

  li1 = tmp[0];
  li0 = tmp[1];

  for (t = 0; t < num_iters - 1; t++) {
    li0[t + 1] = tmp[1][t + 1];
    li0[t] = tmp[0][t];

    for (i = 1; i < t + 1; i++) {
      li1[i] = (li0[i - 1] + li0[i] + li0[i + 1]) / 3.0;
    }

    swap(&li1, &li0);
  }

  t = num_iters - 1;

  for (i = 1; i < t + 1; i++) {
    li1[i] = (li0[i - 1] + li0[i] + li0[i + 1]) / 3.0;
  }

  for (i = 0; i < num_iters; i++) {
    jbi_out[i] = li1[i];
  }
}

static inline void do_topright_hdiam(int num_iters, int pb_size, double **tmp,
  double *jbi_out) {
  int i, l0, t;
  double *li1, *li0;

  l0 = pb_size - num_iters - 2;
  li1 = tmp[0] + l0;
  li0 = tmp[1] + l0;

  for (t = 0; t < num_iters; t++) {
    for (i = pb_size - t - 1; i < pb_size - 1; i++) {
      li1[i-l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
    }
    swap(&li1, &li0);
  }

  for (i = l0 + 1; i < pb_size - 1; i++) {
    jbi_out[i] = li0[i - l0];
  }

}

/*
* task[i][j] is set to 1 if the task on time step i and column j has been
* executed. This function checks if there is no task that has been executed
* without its dependencies being statisfied.
*/
int task_index(uint8_t **tasks, int num_strips, int num_steps) {
  int i, t;

  for (i = 1; i < num_strips; i++) {
    if (tasks[0][i] == 0 && tasks[0][i + 1] == 1) return -1;
  }

  for (t = 1; t < num_steps; t++) {
    for (i = 1; i < num_strips; i++) {
      if ((tasks[i + 1][t - 1] == 0 || tasks[i - 1][t] == 0) &&
          tasks[i][t] == 1)
        return -1;
    }
  }
  return 1;
}

int check_low_iter(int pb_size, int num_stencil_iters) {
  if ((num_stencil_iters >= 16) && (num_stencil_iters <= 64) &&
      (pb_size > num_stencil_iters * (1 << 3))) {
    return 1;
  } else {
    return -1;
  }
}

int check_tilable(int pb_size, int num_stencil_iters) {
  if (pb_size > (T_WIDTH_DBL << 2) && num_stencil_iters > 2 * T_ITERS) {
    return 1;
  } else {
    return -1;
  }
}

int check_default(int pb_size, int num_stencil_iters) {
  if (pb_size > num_stencil_iters && num_stencil_iters > 2) {
    return 1;
  } else {
    return -1;
  }
}
