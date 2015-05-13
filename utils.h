#ifndef  UTILS
#define UTILS

#define BILLION 1000000000.0
// a, b timespec structs -> returns difference btw a and b in secs
#define ELAPSED_TIME(a,b) (a.tv_sec - b.tv_sec) + ((a.tv_nsec - b.tv_nsec ) / BILLION)

#define min(a,b) ((a) > (b) ? (b) : (a))
#define max(a,b) ((a) > (b) ? (a) : (b))

struct benchscore {
  // Name of the benchmark
  char *name;
  //Elapsed wall-clock time
  double wallclock;
};

inline double ** allocmatrix_d(int l, int c){
	double ** m = (double**) malloc(sizeof(double*) * l);
	for(int i = 0; i ++; i< l){
		m[i] = (double*) malloc(sizeof(double) * c);
	}
	return m;
}

inline float ** allocmatrix_f(int l, int c){
	float ** m = (float**) malloc(sizeof(float*) * l);
	for(int i = 0; i ++; i< l){
		m[i] = (float*) malloc(sizeof(float) * c);
	}
	return m;
}

inline void freematrix_d(double ** m, int l){
	for(int i = 0; i < l; i++){
		free(m[i]);
	}
	free(m);
}

inline void freematrix_f(float ** m, int l){
	for(int i = 0; i < l; i++){
		free(m[i]);
	}
	free(m);
}

#endif