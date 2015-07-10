# tiling-strategies

# Jacobi1d

We propose different implementations of the stencil computation for concurrent execution using different shapes of tiles:

0. Naïve implementation : parallelize only the inner loops.
1. Overlapping tiles : the iterations space is partitionned into trapezoïdal tiles that can be computed in parallel. The implementation contains some errors in the calculation a this point.
2. Skewed tiles : paralelogram-shaped tiles, with inter-tile dependences along the horizontal axis and the (-1,1) axis.
3. Half-diamonds, or triangles : for a small number iterations, we need to start the tiles concurrently along the space dimension. We also want to avoid redundant computations, so we use triangles tiles, with two levels of tiles in the computation (see report).

## Compilation & running

'make' and './test -h' for help.

## Viewing half-diamond tiling performance
The aim of this project was to compare the algorithm using a barrier between upper and lower tiles, with static scheduling, and an algorithm relying on tasks. Data from the task implmeentation has been imported in ```data/jacobi1D_starnk_raw.csv```. \\
Run ```stat_from_data.py``` and then ```visualize.py``` in ```data/``` to view some performance graphs. We use ```pandas```, ```numpy```,```matplotlib``` and ```scipy``` packages.
