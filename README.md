# tiling-strategies

# Jacobi1d

Different implementations of one-dimensional jacobi stencil, with different tilings and different approaches :
        - the sequential example
        - "skewed tiles" are paralleogram-shaped tiles adapted for big iteration numbers.
        - "hdiams". For small iterations numbers, we use triangles tiles ( half-diamonds) to enable both a possible concurernt start along the space dimension and parallelism, without redundant computation.
        This tiling is used into classical openmp for loops
                - with a barrier between top and bottom triangles,
                - with grouping tiles, and a barrier into and between groups.
                - with a task-based approach. OpenMp tasks don't seem to perform well here.

## Compilation & running

'make' and './test -h' for help.
