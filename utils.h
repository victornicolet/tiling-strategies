#ifndef  _UTILS_H
#define _UTILS_H

// Cache line size of 64 bytes on most x86
#define  CACHE_LINE_SIZE 64
#define  L1_CACHE_SIZE 6044


#define BILLION 1000000000.0
#define DOUBLE_COMPARISON_THRESHOLD 10e-7
// a, b timespec structs -> returns difference btw a and b in secs
#define ELAPSED_TIME(a,b) (a.tv_sec - b.tv_sec) +                             \
                          ((a.tv_nsec - b.tv_nsec ) / BILLION)

#define min(a,b) ((a) > (b) ? (b) : (a))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define ABS(a) (((a)>0 ) ? (a) : (-a))

struct benchscore {
  // Name of the benchmark
  char *name;
  //Elapsed wall-clock time
  double wallclock;
  int runs;
};

struct benchspec {
  // Name of the benchmark
  char *name;
  // Function to call
  void (*variant)(int, int, double**, struct benchscore *);
  // Function to check alg. arguments
  int (*checkfunc)(int, int);
  // Spatial dimensions
  int dim;
  // Base size
  long size;
  // Iterated stencil multiplier
  int iters;
};

static inline double **
alloc_double_mx(int dim1, int dim2)
{
  int i;
  double ** mx = malloc(dim1 * sizeof ** mx);

  if(mx == NULL){
    fprintf(stderr, "Error while allocating matrix\n");
    return NULL;
  }

  for (i = 0; i < dim1; i++){
    mx[i] = malloc(dim2 * sizeof *mx[i]);

    if(mx[i] == NULL){
      fprintf(stderr, "Error while allocating matrix on line %i\n", i);
      return NULL;
    }
  }

  return mx;
}

static inline void
free_mx(double ** mx, int dim1)
{
  int i;
  for(i = 0; i < dim1; i++){
    free(mx[i]);
  }
  free(mx);
}

int adjust_num(double);

int compare(double *, double *, int);

void swap(void *, void *, size_t);

#endif
