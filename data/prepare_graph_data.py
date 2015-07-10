import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pan
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
from scipy.interpolate import griddata

KB=150


def prep_gdat(write_csv):

        # Data from Starnk
        kpn_raw_dat = np.genfromtxt('jacobi1D_starnk_raw.csv', delimiter=';')

        kpn_rawdf = pan.DataFrame(kpn_raw_dat, columns=['i','s','t'])

        min_i = int(np.min(kpn_rawdf['i']))
        min_s = int(np.min(kpn_rawdf['s']))


        kpn_rawdf['i'] = kpn_rawdf['i'].astype(int)
        kpn_rawdf['s'] = kpn_rawdf['s'].astype(int)
        kpn_rawdf['i'].map(lambda x : int(2<<x))
        kpn_rawdf['s'].map(lambda x : int(2<<x))

        kpn_stats = kpn_rawdf.groupby(['i','s'], sort='yes').agg([np.mean, np.median, np.std])

        print kpn_stats

        kpn_stats.to_csv('jacobi1D_starnk_stats.csv')

        # Data from OpenMP
        omp_stats = np.genfromtxt('jacobi1D_stats.csv', delimiter=',', skip_header=3)
        ompdf = pan.DataFrame(omp_stats, columns=['i','s', 'algo', 'nthreads','omp_mean','omp_median','omp_std'])
        ompdf.drop('nthreads', axis=1, inplace=True)
        ompdf.drop('omp_mean', axis=1, inplace=True)
        ompdf.drop('omp_std', axis=1, inplace=True)
        #Separate data of each algorithm
        omp_hdiam = ompdf.loc[ompdf['algo']==0];
        omp_hdiam_var_tri = ompdf.loc[ompdf['algo']==1];
        omp_hdiam_grp = ompdf.loc[ompdf['algo']==2];
        omp_hdiam_tasks = ompdf.loc[ompdf['algo']==3];

        print omp_hdiam_tasks
        # Merge
        df = pan.merge(omp_hdiam_tasks, kpn_stats, how='inner', on=['i','s'])

        df[['i', 's']] = df[['i', 's']].astype(int)

        #-------------------------------------------------------------------------------
        # Performance : time elapsed for libkpn version / time elapsed for OpenMP version
        df.insert(5, 'spdup_libkpn', 0, allow_duplicates=False)
        df.spdup_libkpn = 100.0 * df.t / df.omp_median;
        df.drop('t', axis=1, inplace=True)


        #-------------------------------------------------------------------------------
        # GFLOPS : for one calculation, we make 2 float add, one division -> 3 FLOPS
        # With doubles, makes 6 FLOPS.

        #nops = 6

        #df.insert(6,'GFLOPS_libkpn', 0, allow_duplicates=False)

        #nops_tot = (df.s * KB - 2) * df.i * nops

        #df.GFLOPS_libkpn = nops_tot / (df.spdup_libkpn * df.seq_t);
        #-------------------------------------------------------------------------------
        if(write_csv):
                df.to_csv('jacobi1D_graph_data.csv', sep=',',index=False)

        return df
