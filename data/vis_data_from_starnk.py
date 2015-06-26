import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pan
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
from scipy.interpolate import griddata

KB=150

raw_dat = np.genfromtxt('jacobi1D_starnk_raw.csv', delimiter=';')

rawdf = pan.DataFrame(raw_dat, columns=['i','s','t'])

max_i = int(np.max(rawdf['i']))
min_i = int(np.min(rawdf['i']))
range_i = int(math.floor(max_i - min_i))
max_s = int(np.max(rawdf['s']))
min_s = int(np.min(rawdf['s']))
range_s = int(math.floor(max_s - min_s))

resulting_data = np.zeros((range_i*range_s, 3))

calculated_means = rawdf.groupby(['i','s'], sort='yes')

for i in range(0, range_i):
    for s in range(0, range_s):
        resulting_data[i * range_s + s][0] = int((2<<(i + min_i)))
        resulting_data[i * range_s + s][1] = int((2<<(s + min_s)))
        resulting_data[i * range_s + s][2] = calculated_means.get_group((i+min_i,s+min_s))['t']._get_numeric_data().mean()

resdf = pan.DataFrame(resulting_data, columns=['i','s','t'])
omp_dat = np.genfromtxt('jacobi1D_viz.csv', delimiter=';', skip_header=1)
ompdf = pan.DataFrame(omp_dat, columns=['i','s','seq_t','spdup_omp','spdup_grps','ignore'])
ompdf.drop('ignore', axis=1, inplace=True)
df = pan.merge(ompdf, resdf, how='inner', on=['i','s'])

df[['i', 's']] = df[['i', 's']].astype(int)

#-------------------------------------------------------------------------------
# Performance : time elapsed for current version / time elapsed for sequential
# version
df.insert(5, 'spdup_libkpn', 0, allow_duplicates=False)
df.spdup_libkpn = 100.0 * df.t / df.seq_t;
df.drop('t', axis=1, inplace=True)


#-------------------------------------------------------------------------------
# GFLOPS : for one calculation, we make 2 float add, one division -> 3 FLOPS
# With doubles, makes 6 FLOPS.

nops = 6

df.insert(6,'GFLOPS_omp', 0, allow_duplicates=False)
df.insert(7,'GFLOPS_grps', 0, allow_duplicates=False)
df.insert(8,'GFLOPS_libkpn', 0, allow_duplicates=False)

nops_tot = (df.s * KB - 2) * df.i * nops

df.GFLOPS_omp = nops_tot / (df.spdup_omp * df.seq_t);
df.GFLOPS_grps = nops_tot / (df.spdup_grps * df.seq_t);
df.GFLOPS_libkpn = nops_tot / (df.spdup_libkpn * df.seq_t);
#-------------------------------------------------------------------------------

print df.columns

df.to_csv('jacobi1D_full.csv', sep=',',index=False)
