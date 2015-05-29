# tiling-strategies

# Jacobi1d

Five implementations of one dimensional Jacobi stencil.
- sequential, and array swapping
- naive OpenMp implementation with parallel for pragmas
- overlapped tiling implementation, tile size is defined by T_WIDTH_DBL_OVERLAP and T_ITERS, a tile fillls a 12kb L1 cache.
- skewed tiling version with task-based approach. The calculation space is tiled with hyperplanes 1:0 along space dimension, 1:1 along time dimension.

## Compilation & running

Jacobi 1D only :
  * $ make jbi1d
  * $ ./jbi1d Nruns Mask [ Width Time iterations ]
  Mask example : 10000 ( overlapped version : naive : skewed : sequential )
  * Some macros defined benchmarks/in jacobi1d.h
