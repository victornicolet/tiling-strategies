#ifndef JACOBI1D
#define JACOBI1D
// App
#define CHECK_ON_SIZE 8
// Cache line size of 64 bytes on most x86
#define  CACHE_LINE_SIZE 64
#define  L1_CACHE_SIZE 6044
// TILE DIMENSIONS -------------------------------
// Iterations within a tile
#define T_ITERS 32
// Fill a cache line with 8 doubles or 16 float
#define T_WIDTH_DBL 32
#define T_WIDTH_FLT 32
// Different values for overlapped version
static int T_WIDTH_DBL_OVERLAP =
 (L1_CACHE_SIZE - T_ITERS * T_ITERS)/(sizeof(double)*T_ITERS);
 // Specifi for diamond tiles
 #define T_WIDTH_DBL_DIAM 8
 #define T_ITERS_DIAM 32
// -------------------------------
#define JBI1D_STENCIL_T(jbi) jbi[t][i] = \
  (jbi[t-1][i] + jbi[t-1][i-1] + jbi[t-1][i+1]) / 3.0
#define JBI1D_STENCIL_SW(jbi) jbi[1][i] = \
  (jbi[0][i-1] + jbi[0][i] + jbi[0][i]) / 3.0
#define JBI1D_STENCIL(lvl1,lvl0) lvl1[i] = \
  (lvl0[i-1] + lvl0[i+1] + lvl0[i]) / 3.0 
//-----------------------------------
#define JBI_INIT(jbi, n) for(j = 0; j < n; j++){\
      jbi[0][j] = cos((double) j );\
      jbi[1][j] = 0;\
    }

void djbi1d_omp_naive_test(int, int, double**, struct benchscore * );
void djbi1d_omp_overlap_test(int, int, double**, struct benchscore * );
void djbi1d_skewed_tiles_test(int, int, double **, struct benchscore * );



#endif /* JACOBI1D */