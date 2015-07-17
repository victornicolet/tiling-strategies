import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pan
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
from scipy.interpolate import griddata



def prep_gdat(write_csv):

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
        #kpn_rawdf['i'] = kpn_rawdf['i'].map(lambda x : int(2<<x))
        kpn_rawdf['s'] = kpn_rawdf['s'].map(lambda x : int(2<<x))

        kpn_naive_rawdf['i'] = kpn_naive_rawdf['i'].astype(int)
        kpn_naive_rawdf['s'] = kpn_naive_rawdf['s'].astype(int)
        #kpn_naive_rawdf['i'] = kpn_naive_rawdf['i'].map(lambda x : int(2<<x))
        kpn_naive_rawdf['s'] = kpn_naive_rawdf['s'].map(lambda x : int(2<<x))

        kpn_merged = pan.merge(kpn_rawdf, kpn_naive_rawdf, how='inner',
                on=['i','s'])

        kpn_stats_t = kpn_merged.groupby(['i','s'],
                sort='yes').agg([np.mean,np.median, np.std])

        kpn_stats = kpn_stats_t.reset_index()

        kpn_stats.columns = [' '.join(col).strip() for col in \
                kpn_stats.columns.values]

        kpn_stats.columns = ['i','s','opt_mean', 'opt_median', 'opt_std',
        'naive_mean', 'naive_median', 'naive_std']


        kpn_stats.to_csv('jacobi1D_starnk_stats.csv')

        # Data from OpenMP
        omp_stats = np.genfromtxt('jacobi1D_stats.csv', delimiter=',',
         skip_header=3)
        ompdf = pan.DataFrame(omp_stats, columns=['i','s', 'algo',
                'nthreads','omp_mean','omp_median','omp_std'])
        ompdf.drop('nthreads', axis=1, inplace=True)
        #Separate data of each algorithm
        omp_hdiam = ompdf.loc[ompdf['algo']==0];
        omp_hdiam_var_tri = ompdf.loc[ompdf['algo']==1];
        omp_hdiam_grp = ompdf.loc[ompdf['algo']==2];
        omp_hdiam_tasks = ompdf.loc[ompdf['algo']==3];

        # Merge
        df = pan.merge(omp_hdiam_var_tri, kpn_stats, how='inner',
                on=['i','s'])

        df.drop('algo', axis=1, inplace=True)

        df[['i', 's']] = df[['i', 's']].astype(int)


        #-------------------------------------------------------------------------------
        # Performance : time elapsed for libkpn version / time elapsed for OpenMP version
        df.insert(2, 'spdup_libkpn_naive_omp', 0, allow_duplicates=False)
        df.spdup_libkpn_naive_omp = df.naive_median / df.omp_median;
        df.insert(2, 'spdup_libkpn_omp', 0, allow_duplicates=False)
        df.spdup_libkpn_omp = df.opt_median / (df.omp_median);
        df.insert(2, 'spdup_libkpn_opt_naive', 0, allow_duplicates=False)
        df.spdup_libkpn_opt_naive = df.opt_median / \
                df.naive_median;

        print df.columns

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
