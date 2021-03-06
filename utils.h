#ifndef _UTILS_H
#define _UTILS_H


/* Number of doubles displayed */
#define DISPLAY_SIZE 8
#ifndef DISPLAY_OFFSET
  #define DISPLAY_OFFSET 1200
#endif
#define DISPLAY_VERBOSE verbose_flag
/*
 * A default number of cores, if nothing works when trying to get the number of
 *  online cores.
 */
#ifndef DEFAULT_PROC_NUM
 #define DEFAULT_PROC_NUM 2
#endif
/* Cache line size of 64 bytes on most x86 */
#define  CACHE_LINE_SIZE 64
#define  L1_CACHE_SIZE 32000

#define KB 125 // Store 125 double precision numbers in a kB
#define BILLION 1000000000.0
#define DOUBLE_COMPARISON_THRESHOLD 10e-7
/* a, b timespec structs -> returns difference btw a and b in secs */
#define ELAPSED_TIME_S(a,b) (a.tv_sec - b.tv_sec) +                             \
                          ((a.tv_nsec - b.tv_nsec ) / BILLION)

#define min(a,b) ((a) > (b) ? (b) : (a))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define ABS(a) (((a)>0 ) ? (a) : (-a))

/* Colors in terminal */
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define KRESET "\033[0m"

#ifdef DEBUG
 #define DEBUG_SIZE 4096
 #define DEBUG_ITER 16
#endif

struct args_dimt
{
  int width, height, iters;
};

struct benchspec
{
  /* Name of the benchmark */
  char *name;
  /*/ Function to call */
  double (*variant)(struct args_dimt, double *, double *);
  /* Function to check alg. arguments */
  int (*checkfunc)(int, int);
  /* Base size */
  long size;
  /* Iterated stencil multiplier */
  int iters;
  /* Spatial dimensions */
  int dim;
};

struct benchspec1d_l
{
  /* Name of the benchmark */
  char *name;
  /*/ Function to call */
  double (*variant)(struct args_dimt, long *, long *);
  /* Function to check alg. arguments */
  int (*checkfunc)(int, int);
  /* Base size */
  long size;
  /* Iterated stencil multiplier */
  int iters;
  /* Spatial dimensions */
  int dim;
};

struct benchspec2d
{
  /* Name of the benchmark */
  char *name;
  /* Function to call */
  double (*variant)(struct args_dimt, double **, double **);
  /* Function to check alg. arguments */
  int (*checkfunc)(int, int, int);
  /* Base size : width and height */
  long width, height;
  /* Iterated stencil multiplier */
  int iters;
};

/* Helpers for allocation of matrixes */

static inline double **
alloc_double_mx(int dim1, int dim2)
{
  int i;
  double ** mx = malloc(dim1 * sizeof ** mx);

  if (mx == NULL) {
    fprintf(stderr, "Error while allocating matrix\n");
    return NULL;
  }

  for (i = 0; i < dim1; i++) {
    mx[i] = malloc(dim2 * sizeof *mx[i]);

    if (mx[i] == NULL) {
      fprintf(stderr, "Error while allocating matrix on line %i\n", i);
      return NULL;
    }
  }

  return mx;
}

static inline long **
alloc_long_mx(int dim1, int dim2)
{
  int i;
  long ** mx = malloc(dim1 * sizeof ** mx);

  if (mx == NULL) {
    fprintf(stderr, "Error while allocating matrix\n");
    return NULL;
  }

  for (i = 0; i < dim1; i++) {
    mx[i] = malloc(dim2 * sizeof *mx[i]);

    if (mx[i] == NULL) {
      fprintf(stderr, "Error while allocating matrix on line %i\n", i);
      return NULL;
    }
  }

  return mx;
}

static inline void
free_mx(void ** mx, int dim1)
{
  int i;
  for (i = 0; i < dim1; i++) {
    free(mx[i]);
  }
  free(mx);
}

static inline double *
alloc_line(int num_elements)
{
  double * line = aligned_alloc(CACHE_LINE_SIZE, num_elements * sizeof(*line));
  if (line == NULL) {
    fprintf(stderr, "Error while allocating line of %i doubles.\n",
            num_elements);
  }
  return line;
}

int adjust_num(double);

long compare(int, double *, double *);

long compare_fast(int, double *, double *);

long compare_l(long *, long *, int);

void find_diffs(int, double *, double *, int *);

void init_data_1d(int, double *);

void init_data_1d_l(int, long *);

void init_data_2d(int, int, double **);

void print_benchspecs(int, struct benchspec *);

void print_check(int, int, double *, double *);

void print_test1d_summary(int, int, double, struct benchspec, double *,
                          double *);

void print_test1d_l_summary(int, int, double, struct benchspec1d_l, long *,
                            long *);

void print_test2d_summary(int, int, double, struct benchspec2d, double **,
                          double **);

void swap(double **, double **);

void why_fopen(int err_no);

#endif
