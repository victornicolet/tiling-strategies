import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pan
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
from scipy.interpolate import griddata


max_gain = 24
max_gain_01=3
watch_size = [32768, 65536, 131072]

# Data from Libkpn - optimized version
kpn_raw_dat = np.genfromtxt('jacobi1D_starnk_raw.csv', delimiter=';')

# Data from Libkpn - naive version
kpn_raw_naive = np.genfromtxt('jacobi1D_starnk_naive_raw.csv',
        delimiter=';')

kpn_rawdf = pan.DataFrame(kpn_raw_dat, columns=['i','s','t_opt'])
kpn_naive_rawdf = pan.DataFrame(kpn_raw_naive,
        columns=['i','s','t_naive'])

min_i = int(np.min(kpn_rawdf['i']))
min_s = int(np.min(kpn_rawdf['s']))
min_n_i = int(np.min(kpn_naive_rawdf['i']))
min_n_s = int(np.min(kpn_naive_rawdf['s']))


kpn_rawdf['i'] = kpn_rawdf['i'].astype(int)
kpn_rawdf['s'] = kpn_rawdf['s'].astype(int)

kpn_rawdf['s'] = kpn_rawdf['s'].map(lambda x : int(2<<x))

kpn_naive_rawdf['i'] = kpn_naive_rawdf['i'].astype(int)
kpn_naive_rawdf['s'] = kpn_naive_rawdf['s'].astype(int)

kpn_naive_rawdf['s'] = kpn_naive_rawdf['s'].map(lambda x : int(2<<x))

print kpn_naive_rawdf
print kpn_rawdf

kpn_merged = pan.merge(kpn_rawdf, kpn_naive_rawdf, how='inner',
        on=['i','s'])

kpn_stats_t = kpn_merged.groupby(['i','s'],
        sort='yes').agg([np.mean,np.median, np.std])


kpn_stats = kpn_stats_t.reset_index()

kpn_stats.columns = [' '.join(col).strip() for col in \
        kpn_stats.columns.values]

kpn_stats.columns = ['i','s','opt_mean', 'opt_median', 'opt_std',
'naive_mean', 'naive_median', 'naive_std']


kpn_stats.to_csv('jacobi1D_kpnonly_stats_.csv')


#-------------------------------------------------------------------------------
# Performance : time elapsed for libkpn version / time elapsed for OpenMP version
kpn_stats.insert(2, 'spdup_libkpn_opt_naive', 0, allow_duplicates=False)
kpn_stats.spdup_libkpn_opt_naive = kpn_stats.opt_median / \
        kpn_stats.naive_median;

my_data = kpn_stats.as_matrix()

# note this: you can skip rows!
X = my_data[:,0]
Y = np.vectorize(math.log)(my_data[:,1], 2)
Z1 = np.vectorize(lambda x: 1/x)(my_data[:,2 ])

index = np.where(Y == math.log(watch_size[0], 2))
index2 = np.where(Y == math.log(watch_size[1], 2))
index3 = np.where(Y == math.log(watch_size[2], 2))

Xb = X[index]

Zb1 = Z1[index]
Zb1_2 = Z1[index2]
Zb1_3 = Z1[index3]

# comparison for specific problem size
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)

ax2.set_xlabel('Iterations')
ax2.set_ylabel('Speedup opt / naive')

l1, = ax2.plot(Xb,Zb1, color='red')
l2, = ax2.plot(Xb,Zb1_2, color='blue')
l3, = ax2.plot(Xb,Zb1_3, color='green')

fig2.legend([l1,l2,l3],['32768 Kb','65536 Kb', '131072 Kb'],
        loc='bottom', ncol=3)


plt.show()