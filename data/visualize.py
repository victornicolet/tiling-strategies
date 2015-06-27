import csv
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
from scipy.interpolate import griddata

max_gain = 24
watch_size = 32768

fig = plt.figure(figsize=plt.figaspect(0.5))
fig.text(1, 1, "Benchmark : jacobi1d on 12 cores", fontdict=None)
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
ax3 = fig.add_subplot(2, 2, 3, projection='3d')
ax4 = fig.add_subplot(2, 2, 4)

ax1.set_xlabel('Iterations');
ax1.set_ylabel('log(size_kB)');
ax1.set_zlabel('Speedup - not grouped');

ax2.set_xlabel('Iterations');
ax2.set_ylabel('log(size_kB)');
ax2.set_zlabel('Speedup - grouped');

ax3.set_xlabel('Iterations');
ax3.set_ylabel('log(size_kB)');
ax3.set_zlabel('Speedup - starnk');

ax4.set_xlabel('Iterations, size : 32768 kB')
ax4.set_ylabel('Speedup - not grouped vs. grouped');
# note this: you can skip rows!
my_data = np.genfromtxt('jacobi1D_full.csv', delimiter=',',skiprows=1)
X = my_data[:,0]
Y = np.vectorize(math.log)(my_data[:,1], 2)
Z1 = np.vectorize(lambda x: 100/x)(my_data[:,3])
Z2 = np.vectorize(lambda x: 100/x)(my_data[:,4])
Z3 = np.vectorize(lambda x: 100/x)(my_data[:,5])

index = np.where(Y == math.log(watch_size, 2))
print index
Xb = X[index]
Zb1 = Z1[index]
Zb2 = Z2[index]
Zb3 = Z3[index]

xi = np.linspace(X.min(),X.max(),100)
yi = np.linspace(Y.min(),Y.max(),100)
# VERY IMPORTANT, to tell matplotlib how is your data organized
zi1 = griddata((X, Y), Z1, (xi[None,:], yi[:,None]), method='cubic')
zi2 = griddata((X, Y), Z2, (xi[None,:], yi[:,None]), method='cubic')
zi3 = griddata((X, Y), Z3, (xi[None,:], yi[:,None]), method='cubic')

xig, yig = np.meshgrid(xi, yi)

ax1.set_zlim([0, max_gain]);
ax2.set_zlim([0, max_gain]);
ax3.set_zlim([0, max_gain]);

surf1 = ax1.plot_surface(xig, yig, zi1,
        linewidth=0, color='b')
surf2 = ax2.plot_surface(xig, yig, zi2,
        linewidth=0, color='r')
surf3 = ax3.plot_surface(xig, yig, zi3,
        linewidth=0, color='g')

ax4.set_ylim([0, max_gain]);
ax4.plot(Xb, Zb1, 'ro', color='b')
ax4.plot(Xb, Zb2, 'ro', color='r')
ax4.plot(Xb, Zb3, 'ro', color='g')
ax4.axhline(1, color='r')

fig2 = plt.figure(figsize=plt.figaspect(0.5))
ax = fig2.add_subplot(1,1,1);

omp_data = np.genfromtxt('jacobi1D_viz.csv', delimiter=';', skiprows=1)
Xomp = omp_data[:,0]
Yomp = np.vectorize(math.log)(my_data[:,1], 2)
Zomp1 = omp_data[:,3]
Zomp2 = omp_data[:,4]
indx = np.where(Yomp == math.log(watch_size, 2))
Xompb = Xomp[index]
Zomp1b = Zomp1[index]
Zomp2b= Zomp2[index]

ax.plot(Xb, Zomp1b, 'ro', color='r')
ax.plot(Xb, Zomp2b, 'ro', color='b')

plt.show()
