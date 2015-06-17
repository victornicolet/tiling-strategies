import csv
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
from scipy.interpolate import griddata

fig = plt.figure(figsize=plt.figaspect(0.5))
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
ax3 = fig.add_subplot(2, 2, 3)

ax1.set_xlabel('Iterations');
ax1.set_ylabel('log(size_kB)');
ax1.set_zlabel('Speedup - not grouped');

ax2.set_xlabel('Iterations');
ax2.set_ylabel('log(size_kB)');
ax2.set_zlabel('Speedup - grouped');

ax3.set_xlabel('Iterations')
ax3.set_ylabel('Speedup - not grouped vs. grouped');
# note this: you can skip rows!
my_data = np.genfromtxt('data/jacobi1D_viz.csv', delimiter=';',skiprows=1)
X = my_data[:,0]
Y = np.vectorize(math.log)(my_data[:,1], 2)
Z1 = np.vectorize(lambda x: 100/x)(my_data[:,3])
Z2 = np.vectorize(lambda x: 100/x)(my_data[:,4])

index = np.where(Y == math.log(32768, 2))
print index
Xb = X[index]
Zb1 = Z1[index]
Zb2 = Z2[index]

xi = np.linspace(X.min(),X.max(),100)
yi = np.linspace(Y.min(),Y.max(),100)
# VERY IMPORTANT, to tell matplotlib how is your data organized
zi1 = griddata((X, Y), Z1, (xi[None,:], yi[:,None]), method='cubic')
zi2 = griddata((X, Y), Z2, (xi[None,:], yi[:,None]), method='cubic')

xig, yig = np.meshgrid(xi, yi)

surf1 = ax1.plot_surface(xig, yig, zi1,
        linewidth=0)
surf2 = ax2.plot_surface(xig, yig, zi2,
        linewidth=0, color='r')
ax3.set_ylim([0,16]);
ax3.plot(Xb, Zb1, 'ro')
ax3.plot(Xb, Zb2, 'ro', color='b')

plt.show()
