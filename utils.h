#ifndef  UTILS
#define UTILS

#define BILLION 1000000000.0
// a, b timespec structs -> returns difference btw a and b in secs
#define ELAPSED_TIME(a,b) (a.tv_sec - b.tv_sec) + ((a.tv_nsec - b.tv_nsec ) / BILLION)

#define min(a,b) ((a) > (b) ? (b) : (a))
#define max(a,b) ((a) > (b) ? (a) : (b))

#define ALLOC_MX(m, type, dim1, dim2) type ** (m) = \
	(type **) malloc(sizeof(type *) * (dim1)); \
	for(int i = 0; i < (dim1); i++){\
	m[i] = (type *) malloc(sizeof(type) * dim2);\
	}
	
#define FREE_MX(m, dim1) for(int i = 0; i < (dim1); i++){free(m[i]);}free(m);

struct benchscore {
  // Name of the benchmark
  char *name;
  //Elapsed wall-clock time
  double wallclock;
};

#endif