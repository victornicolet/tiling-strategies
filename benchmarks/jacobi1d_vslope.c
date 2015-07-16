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

#include <omp.h>

static struct timespec tend;
static struct timespec tbegin;

/* Functions describing different tasks in the computation
  Always inlined in the main body */
static inline void do_top_hdiam(int, int, int, int, double **, double *)
    __attribute__((always_inline));
static inline void do_base_hdiam(int, int, int, int,  double *, double **)
    __attribute__((always_inline));
static inline void do_topleft_hdiam(int, int, double **, double *)
    __attribute__((always_inline));
static inline void do_topright_hdiam(int, int, int, double **, double *)
    __attribute__((always_inline));


int get_slope(int num_iters) {
  return max(
    ((int) floor(L1_CACHE_SIZE / (num_iters * sizeof(double) * 4))) / 4 * 4,
     1);
}


static void djbi1d_hdiam_vslope(int pb_size, int num_iters, double *jbi_in,
  double *jbi_out, double **tmp) {
  int num_tiles, slope, tile_no, tile_size;
  /* Slope caclulation : fit the tile in L1 cache */
  slope = get_slope(num_iters);
  tile_size = 2 * slope * num_iters;
  num_tiles = pb_size / tile_size;

    /* First loop : base down tiles */
#pragma omp parallel for schedule(static) shared(tmp) private(tile_no)
  for (tile_no = 0; tile_no < num_tiles; tile_no++) {
    do_base_hdiam(num_iters, pb_size, tile_no, slope, jbi_in, tmp);
  }

  /* Half-tile at beggining, in front of parallel loop */
 do_topleft_hdiam(num_iters, slope, tmp, jbi_out);

  /* Second loop : tip down tiles */
#pragma omp parallel for schedule(static) shared(tmp) private(tile_no)
  for (tile_no = 1; tile_no < num_tiles; tile_no++) {
    #ifdef TEST_LOCALITY
        tiles[1][tile_no] = omp_get_thread_num();
    #endif
    do_top_hdiam(num_iters, pb_size, tile_no, slope, tmp, jbi_out);
  }

  do_topright_hdiam(num_iters, pb_size, slope, tmp, jbi_out);
  free_mx((void **)tmp, 2);
}


double djbi1d_hdiam_vslope_t(struct args_dimt args, double *jbi_in,
  double *jbi_out) {

  double **tmp;

  /* Inter-tiles buffer, twice the problem size */
  tmp = aligned_alloc(CACHE_LINE_SIZE, 2 * sizeof(*tmp));

  if (tmp == NULL) {
    printf("Allocation of tmp failed !  ... Aborting.\n");
    return -1.0;
  }

  tmp[0] = aligned_alloc(CACHE_LINE_SIZE, args.width * sizeof(*(tmp[0])));
  tmp[1] = aligned_alloc(CACHE_LINE_SIZE, args.width * sizeof(*(tmp[1])));

  if (tmp[0] == NULL || tmp[1] == NULL) {
    printf("Allocation of tmp lines failed !  ... Aborting\n");
    free(tmp);
    return -1.0;
  }

  clock_gettime(CLOCK_MONOTONIC, &tbegin);
  djbi1d_hdiam_vslope(args.width, args.iters, jbi_in, jbi_out, tmp);
  clock_gettime(CLOCK_MONOTONIC, &tend);

  free(tmp[0]);
  free(tmp[1]);
  free(tmp);

  return ELAPSED_TIME(tend, tbegin);
}

static inline void do_base_hdiam(int num_iters, int pb_size, int tile_no,
 int slope, double *jbi_in, double **tmp) {
  int i, t;
  int r0, l0, r, l;
  int tile_base_sz = 2 * num_iters * slope;
  double *li1, *li0;


  /* Initial values */
  l0 = max(tile_no * tile_base_sz, 0);
  r0 = min(l0 + tile_base_sz, pb_size - 1);

  li1 = tmp[1] + l0;
  li0 = tmp[0] + l0;

  for (i = l0; i < r0; i++) {
    li0[i - l0] = jbi_in[i];
    li1[i - l0] = 0.0;
  }

  for (t = 0; t < num_iters - 1; t++) {
    l = max(l0 + slope * (t + 1), 1);
    r = min(l0 + tile_base_sz - slope * (t + 1), pb_size);

    for (i = l; i < r; i++) {
      li1[i - l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
    }

    swap(&li0, &li1);
  }
}

static inline void do_top_hdiam(int num_iters, int pb_size, int tile_no,
  int slope, double **tmp, double *jbi_out) {
  int t, i;
  int x0, l0, r0, l, r;
  int tile_base_sz = 2 * slope * num_iters;

  x0 = tile_no * tile_base_sz;
  l0 = max(x0 - slope * num_iters - 1, 0);
  r0 = min(x0 + slope * num_iters + 1, pb_size - 1);

  double *li1, *li0;

  li1 = tmp[1] + l0;
  li0 = tmp[0] + l0;

  for (t = 0; t < num_iters; t++) {

    l = max(x0 - slope * (t + 1), 1);
    r = min(x0 + slope * (t + 1), pb_size - 2);

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

static inline void do_topleft_hdiam(int num_iters, int slope, double **tmp,
  double *jbi_out) {
  int i, t;
  double *li1, *li0;

  li1 = tmp[1];
  li0 = tmp[0];

  for (t = 0; t < num_iters - 1; t++) {
    for (i = 1; i < slope * (t + 1); i++) {
      li1[i] = (li0[i - 1] + li0[i] + li0[i + 1]) / 3.0;
    }
    swap(&li1, &li0);
  }

  for (i = 1; i < slope * (t + 1); i++) {
    li1[i] = (li0[i - 1] + li0[i] + li0[i + 1]) / 3.0;
  }

  for (i = 0; i < slope * num_iters; i++) {
    jbi_out[i] = li1[i];
  }
}

static inline void do_topright_hdiam(int num_iters, int pb_size, int slope,
  double **tmp, double *jbi_out) {
  int i, l0, t;
  double *li1, *li0;

  l0 = pb_size - slope * num_iters - 2;
  li1 = tmp[1] + l0;
  li0 = tmp[0] + l0;

  for (t = 0; t < num_iters; t++) {
    for (i = pb_size - slope * (t + 1); i < pb_size - 1; i++) {
      li1[i-l0] = (li0[i - 1 - l0] + li0[i - l0] + li0[i + 1 - l0]) / 3.0;
    }
    swap(&li1, &li0);
  }

  for (i = l0 + 1; i < pb_size - 1; i++) {
    jbi_out[i] = li0[i - l0];
  }

}
