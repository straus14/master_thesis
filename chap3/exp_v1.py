from __future__ import division
import numpy as np
import mpmath as mp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save


fig = plt.figure()
ax = Axes3D(fig)

#x = np.arange(0,15,5)
#y = np.arange(0,15,5)
outfile = np.load('rate_coverage_two_community.npz')
x = outfile['threshold_user1']/10**6
y = outfile['threshold_user2']/10**6
x,y = np.meshgrid(x,y)

z = outfile['rate_cov_comm']

surf = ax.plot_surface(x, y, z, rstride= 1, cstride= 1, cmap='winter')
fig.colorbar(surf,shrink = 0.25, aspect = 10)
ax.set_xlabel('$\\rho^{1} (Mbps)$')
ax.set_ylabel('$\\rho^{2} (Mbps)$')
ax.set_zlabel('$\mathcal{R}$ - Rate Coverage')
plt.title("$\lambda_{u1}=25$ $\lambda_{u2}=25$")
plt.show()
#plt.savefig('fig1.eps')
#tikz_save('fig1.tex')


########################################################################################################################
          # example to save the file
########################################################################################################################
'''
th_u1 = np.arange(0,4,0.2)
th_u2 = np.arange(0,4,0.2)
k = np.zeros((3,3))
np.savez('test',k=k,th_u1=th_u1,th_u2=th_u2)

h = np.load('test.npz')
print(h['k'])
print(h['th_u1'])
print(h['th_u2'])'''



########################################################################################################################
          # slicing and appending to the array
########################################################################################################################
'''th_u1 = np.arange(0,4,0.2)
th_u2 = np.arange(0,4,0.2)
h = th_u1.size
print(h)
k = np.zeros((5,5))

a1 = th_u1[1:]
print(th_u1)
print(a1)
k1 = k[1:,1:]
print(k1)'''