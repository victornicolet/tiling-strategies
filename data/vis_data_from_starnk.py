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
ompdf = pan.DataFrame(omp_dat, columns=['i','s','seq_t','spdup','spdup_groups','ignore'])
ompdf.drop('ignore', axis=1, inplace=True)
df = pan.merge(ompdf, resdf, how='inner', on=['i','s'])

df[['i', 's']] = df[['i', 's']].astype(int)


df.t = 100.0 * df.t / df.seq_t;
print df

df.to_csv('jacobi1D_full.csv', sep=',',index=False)
