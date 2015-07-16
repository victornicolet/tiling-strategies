import csv
import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

from numpy import genfromtxt
from scipy.interpolate import griddata

import prepare_graph_data as pgd


my_data = pgd.prep_gdat(True).as_matrix()

max_gain = 24
max_gain_01=3
watch_size = [32768, 65536, 131072]



fig = plt.figure(figsize=plt.figaspect(0.5))
fig.text(1, 1, "Benchmark : jacobi1d on 12 cores", fontdict=None)
ax1 = fig.add_subplot(1, 1, 1, projection='3d')

ax1.set_xlabel('Iterations');
ax1.set_ylabel('log(size_kB)');
ax1.set_zlabel('Speedup - opt / naive');

# note this: you can skip rows!
X = my_data[:,0]
Y = np.vectorize(math.log)(my_data[:,1], 2)
Z1 = np.vectorize(lambda x: 1/x)(my_data[:,2 ])
Zomp = np.vectorize(lambda x: 1/x)(my_data[:,3 ])
Zomp_n = np.vectorize(lambda x: 1/x)(my_data[:,4 ])

index = np.where(Y == math.log(watch_size[0], 2))
index2 = np.where(Y == math.log(watch_size[1], 2))
index3 = np.where(Y == math.log(watch_size[2], 2))

Xb = X[index]

Zb1 = Z1[index]
Zb1_2 = Z1[index2]
Zb1_3 = Z1[index3]

Zomp1 = Zomp[index]
Zomp2 = Zomp[index2]
Zomp3 = Zomp[index3]

Zomp_n1 = Zomp_n[index]
Zomp_n2 = Zomp_n[index2]
Zomp_n3 = Zomp_n[index3]

xi = np.linspace(X.min(),X.max(),100)
yi = np.linspace(Y.min(),Y.max(),100)
# VERY IMPORTANT, to tell matplotlib how is your data organized
zi1 = griddata((X, Y), Z1, (xi[None,:], yi[:,None]), method='cubic')

xig, yig = np.meshgrid(xi, yi)

# ax1 : plot 3d surface of comparison between libkpn tasks with two
# tiles / one task per tile

ax1.set_zlim([0, max_gain_01]);

surf1 = ax1.plot_surface(xig, yig, zi1,
        linewidth=0, color='b')

point  = np.array([0, 0, 1])
normal = np.array([0, 0, 1])
d = -point.dot(normal)
xx, yy = np.meshgrid(np.arange(X.min(), X.max()), np.arange(Y.min(),
 Y.max()))
z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]
plane = ax1.plot_surface(xx,yy,z,linewidth=0, color='r')

# comparison for specific problem size
fig2 = plt.figure()
ax2 = fig2.add_subplot(3,1,1)
ax3 = fig2.add_subplot(3,1,2)
ax4 = fig2.add_subplot(3,1,3)

ax2.set_xlabel('Iterations')
ax2.set_ylabel('Speedup opt/separate tiles')

ax3.set_xlabel('Iterations')
ax3.set_ylabel('Speedup Libkpn opt / OpenMP')

ax4.set_xlabel('Iterations')
ax4.set_ylabel('Speedup Libkpn naive / OpenMP')

l1, = ax2.plot(Xb,Zb1, color='red')
l2, = ax2.plot(Xb,Zb1_2, color='blue')
l3, = ax2.plot(Xb,Zb1_3, color='green')

ax3.plot(Xb,Zomp1, color='red')
ax3.plot(Xb,Zomp2, color='blue')
ax3.plot(Xb,Zomp3, color='green')

ax4.plot(Xb,Zomp_n1, color='red')
ax4.plot(Xb,Zomp_n2, color='blue')
ax4.plot(Xb,Zomp_n3, color='green')

fig2.legend([l1,l2,l3],['32768 Kb','65536 Kb', '131072 Kb'],
        loc='bottom', ncol=3)


plt.show()
