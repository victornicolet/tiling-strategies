import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pan
from numpy import genfromtxt

omp_dat = np.genfromtxt('jacobi1D_raw.csv', delimiter=',', skip_header=1)

ompdf = pan.DataFrame(omp_dat, columns=['i','s','algo','threads','time'])
ompdf.convert_objects(convert_numeric=True)
ompdf[['i','s','algo','threads']] = ompdf[['i','s','algo','threads']].astype(int)

gmes = ompdf.groupby(['i','s','algo','threads'], sort='yes')

stats = gmes.agg([np.mean, np.median, np.std])

stats.convert_objects(convert_numeric=True)

stats.to_csv("jacobi1D_stats.csv", sep =',')