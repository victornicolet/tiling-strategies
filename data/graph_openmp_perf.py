import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pan
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
from scipy.interpolate import griddata


watch_size = [32768, 65536, 131072]

omp_dat = np.genfromtxt('jacobi1D_stats.csv', skiprows=3, delimiter =',')

I = omp_dat[:, 0]
S = omp_dat[:, 1]
T_median = omp_dat[:, 5]

index = np.where(S == watch_size[0])
index2 = np.where(S == watch_size[1])
index3 = np.where(S == watch_size[2])

I_1 = I[index]
I_3 = I[index2]
I_2 = I[index3]

T_median_1 = T_median[index]
T_median_3 = T_median[index2]
T_median_2 = T_median[index3]

li = T_median_1.__len__()

# T_median_1 = T_median_1 / T_median_1[-1]
# T_median_2 = T_median_2 / T_median_2[-1]
# T_median_3 = T_median_3 / T_median_3[-1]


l1, = plt.plot(I_1, T_median_1)
l2, = plt.plot(I_2, T_median_2)
l3, = plt.plot(I_3, T_median_3)

plt.legend([l1,l2,l3],['32768 Kb','65536 Kb', '131072 Kb'],
        loc='bottom', ncol=3)

plt.show()
