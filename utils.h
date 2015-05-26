#ifndef  UTILS
#define UTILS

// Cache line size of 64 bytes on most x86
#define  CACHE_LINE_SIZE 64
#define  L1_CACHE_SIZE 6044


#define BILLION 1000000000.0
// a, b timespec structs -> returns difference btw a and b in secs
#define ELAPSED_TIME(a,b) (a.tv_sec - b.tv_sec) + ((a.tv_nsec - b.tv_nsec ) / BILLION)

#define min(a,b) ((a) > (b) ? (b) : (a))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define ABS(a) (((a)>0 ) ? (a) : (-a))
#define ALLOC_MX(m, type, dim1, dim2) type ** m = \
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
};


int adjust_num(double num) {
    double low_bound = 1e7;
    double high_bound = low_bound*10;
    double adjusted = num;
    int is_negative = (num < 0);
    if(num == 0) {
        return 0;
    }
    if(is_negative) {
        adjusted *= -1;
    }
    while(adjusted < low_bound) {
        adjusted *= 10;
    }
    while(adjusted >= high_bound) {
        adjusted /= 10;
    }
    if(is_negative) {
        adjusted *= -1;
    }
    return round(adjusted);
}

int compare(double * t1, double * t2, int n){
  for(int i = 0; i < n; i++){
    if(adjust_num(t1[i]) != adjust_num(t2[i])){
      return 0;
    }
  }
  return 1;
}

void swap(void *a, void *b, size_t size) {
  enum { threshold = (1 << 7) };
  if (size <= threshold) {
    char temp[threshold];

    memcpy(temp, b,    size);
    memcpy(b,    a,    size);
    memcpy(a,    temp, size);
  }
  else {
    void* temp = aligned_alloc(CACHE_LINE_SIZE, size);

    memcpy(temp, b,    size);
    memcpy(b,    a,    size);
    memcpy(a,    temp, size);

    free(temp);
  }
}

#endif