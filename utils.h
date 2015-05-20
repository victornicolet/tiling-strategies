#ifndef  UTILS
#define UTILS

#define BILLION 1000000000.0
// a, b timespec structs -> returns difference btw a and b in secs
#define ELAPSED_TIME(a,b) (a.tv_sec - b.tv_sec) + ((a.tv_nsec - b.tv_nsec ) / BILLION)

#define min(a,b) ((a) > (b) ? (b) : (a))
#define max(a,b) ((a) > (b) ? (a) : (b))

#define ALLOC_MX(m, type, dim1, dim2) type ** (m) = \
	(type **) malloc(sizeof(type *) * (dim1)); \
  if(m == NULL){\
    fprintf(stderr, "ALLOC_MX:Error while allocating 2D array\n");\
  }\
	for(int i = 0; i < (dim1); i++){\
	m[i] = (type *) malloc(sizeof(type) * dim2);\
  if(m[i] == NULL){\
  fprintf(stderr, "ALLOC_MX:Error while allocating 2D array at line %i\n", i); \
  }\
	}
	
#define FREE_MX(m, dim1) if( m == NULL){ \
      fprintf(stderr, "FREE_MX:Error freeing NULL ! \n"); \
    } else { \
      for(int i = 0; i < (dim1); i++){ \
        if(m[i] == NULL){ \
          fprintf(stderr, "FREE_MX:Error freeing NULL at line %i\n", i);\
        } else {\
          free(m[i]);\
        }\
      }\
      free(m);\
    }

#define SWAP(l1, l2, tmp)  { tmp = l1; l1 = l2; l2 = tmp; }
struct benchscore {
  // Name of the benchmark
  char *name;
  //Elapsed wall-clock time
  double wallclock;
};


#endif