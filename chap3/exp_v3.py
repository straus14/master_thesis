__author__ = 'straus14'
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
#-----------------------------------------------------------------------------------------------------------------------
          #figures for chapter 3
#-----------------------------------------------------------------------------------------------------------------------
#                     0       1            2              3        4          5          6
#                   [M-rat, K-tier, openaccess_flag, Pow(dBm), density,     alpha,    Bias(dB)]
m_1_k_1 =  np.array([1,       1,           True,          53,        1,         4,         0],float)
m_1_k_2 =  np.array([1,       2,           True,          33,        5,         4,         0],float)
m_1_k_3 = np.array([1,        3,           True,          23,       10,         4,         0],float)
mk_mat1 = np.array([ m_1_k_1, m_1_k_2, m_1_k_3])

#fig = plt.figure()
#ax = Axes3D(fig)

#--------------------------------------fig1a----------------------------------------------------------------
'''
ref_file = np.load('fig1a.npz')
x = ref_file['rate_th1']/10**6
y = ref_file['rate_th2']/10**6
x,y = np.meshgrid(x,y)

z = ref_file['rate_cov_result_fig1a']

surf = ax.plot_surface(x, y, z, rstride= 1, cstride= 1, cmap='rainbow')
fig.colorbar(surf,shrink = 0.25, aspect = 10)
ax.set_xlabel('$\\rho_{1} (\\times 10^{6})$')
ax.set_ylabel('$\\rho_{2} (\\times 10^{6})$')
ax.set_zlabel('$\mathcal{R}_{\kappa}$ - Rate Coverage')
plt.title("$\lambda_{u_1}=25$ $\lambda_{u_2}=25$")
plt.show()
'''
'''
#--------------------------------------fig1b----------------------------------------------------------------
ref_file = np.load('fig1b.npz')
x = ref_file['rate_th1']/10**6
y = ref_file['rate_th2']/10**6
x,y = np.meshgrid(x,y)

z = ref_file['rate_cov_result_fig1b']

surf = ax.plot_surface(x, y, z, rstride= 1, cstride= 1, cmap='rainbow')
fig.colorbar(surf,shrink = 0.25, aspect = 10)
ax.set_xlabel('$\\rho_{2} (\\times 10^{6})$')
ax.set_ylabel('$\\rho_{1} (\\times 10^{6})$')
ax.set_zlabel('$\mathcal{R}_{\kappa}$ - Rate Coverage')
plt.title("$\lambda_{u_1}=38$ $\lambda_{u_2}=12$")
plt.show()
#--------------------------------------fig1b----------------------------------------------------------------
'''
'''
#--------------------------------------fig2a----------------------------------------------------------------
ref_file = np.load('fig2a.npz')
x = ref_file['lamda_u1_ary']
y = ref_file['lamda_u2_ary']
x,y = np.meshgrid(x,y)

z = ref_file['rate_cov_result_fig2a']

surf = ax.plot_surface(x, y, z, rstride= 1, cstride= 1, cmap='rainbow')
fig.colorbar(surf,shrink = 0.25, aspect = 10)
ax.set_xlabel('$\lambda_{u_2}$')
ax.set_ylabel('$\lambda_{u_1}$')
ax.set_zlabel('$\mathcal{R}_{\kappa}$ - Rate Coverage')
plt.title("$\\rho_{1}=0.5\\times10^{6}$ $\\rho_{2}=2\\times10^{6}$")
plt.show()
#--------------------------------------fig2a----------------------------------------------------------------
'''
'''
#--------------------------------------fig2b----------------------------------------------------------------
ref_file = np.load('fig2b.npz')
x = ref_file['lamda_u1_ary']
y = ref_file['lamda_u2_ary']
x,y = np.meshgrid(x,y)

z = ref_file['rate_cov_result_fig2b']

surf = ax.plot_surface(x, y, z, rstride= 1, cstride= 1, cmap='rainbow')
fig.colorbar(surf,shrink = 0.25, aspect = 10)
ax.set_xlabel('$\lambda_{u_2}$')
ax.set_ylabel('$\lambda_{u_1}$')
ax.set_zlabel('$\mathcal{R}_{\kappa}$ - Rate Coverage')
plt.title("$\\rho_{1}=0.5\\times10^{6}$ $\\rho_{2}=4\\times10^{6}$")
plt.show()
#--------------------------------------fig2a----------------------------------------------------------------
'''
'''
#--------------------------------------fig3a----------------------------------------------------------------
ref_file1 = np.load('fig1a.npz')
ref_file2 = np.load('fig3a.npz')
x = ref_file1['rate_th1']/10**6
y = ref_file1['rate_th2']/10**6
x,y = np.meshgrid(x,y)

z = abs(ref_file1['rate_cov_result_fig1a']-ref_file2['rate_cov_result_fig3a'])*100/(ref_file2['rate_cov_result_fig3a'])

surf = ax.plot_surface(x, y, z, rstride= 1, cstride= 1, cmap='rainbow')
fig.colorbar(surf,shrink = 0.25, aspect = 10)
ax.set_xlabel('$\\rho_{2} (\\times 10^{6})$')
ax.set_ylabel('$\\rho_{1} (\\times 10^{6})$')
ax.set_zlabel('G - Gain')
plt.title("$\lambda_{u_1}=25$ $\lambda_{u_2}=25$")
plt.show()
#--------------------------------------fig3a----------------------------------------------------------------
'''
'''
#--------------------------------------fig3b----------------------------------------------------------------
ref_file1 = np.load('fig1b.npz')
ref_file2 = np.load('fig3b.npz')
x = ref_file1['rate_th1']/10**6
y = ref_file1['rate_th2']/10**6
x,y = np.meshgrid(x,y)

z = abs(ref_file1['rate_cov_result_fig1b']-ref_file2['rate_cov_result_fig3b'])*100/(ref_file2['rate_cov_result_fig3b'])

surf = ax.plot_surface(x, y, z, rstride= 1, cstride= 1, cmap='rainbow')
fig.colorbar(surf,shrink = 0.25, aspect = 10)
ax.set_xlabel('$\\rho_{2} (\\times 10^{6})$')
ax.set_ylabel('$\\rho_{1} (\\times 10^{6})$')
ax.set_zlabel('G - Gain')
plt.title("$\lambda_{u_1}=38$ $\lambda_{u_2}=12$")
plt.show()
#--------------------------------------fig3b----------------------------------------------------------------
'''
'''
#--------------------------------------fig5a----------------------------------------------------------------
ref_file = np.load('fig5a.npz')
kappa_array = ref_file['kappa_array']
rate_cov1 = ref_file['rate_cov1']
rate_cov2 = ref_file['rate_cov2']
rate_cov = ref_file['rate_cov']
plt.plot(kappa_array,rate_cov1,'b-*',label='$\mathcal{R}_{1}$')
plt.plot(kappa_array,rate_cov2,'r-*',label='$\mathcal{R}_{2}$')
plt.plot(kappa_array,rate_cov,'g-*',label='$\mathcal{R}_{\kappa}$')
plt.xlabel('$\kappa$')
plt.ylabel('rate coverage')
plt.legend(loc='best')
plt.show()
#--------------------------------------fig5a----------------------------------------------------------------
'''
'''
#--------------------------------------fig6a----------------------------------------------------------------
ref_file = np.load('fig6a.npz')
plt.plot(ref_file['lamda_sc'],ref_file['rate_cov'][0],'b-*',label='$\lambda_{1} = 0$')
plt.plot(ref_file['lamda_sc'],ref_file['rate_cov'][1],'g-*',label='$\lambda_{1} = 5$')
plt.plot(ref_file['lamda_sc'],ref_file['rate_cov'][2],'r-*',label='$\lambda_{1} = 10$')
plt.legend(loc='best')
plt.xlabel('$\lambda_{3}$')
plt.ylabel('rate coverage $\mathcal{R}_{\kappa_p}$')
plt.show()
#--------------------------------------fig6a----------------------------------------------------------------
#print(ref_file.files)
'''
'''
#--------------------------------------fig7a----------------------------------------------------------------
ref_file = np.load('fig7a.npz')
plt.plot(ref_file['lamda_sc'],ref_file['rate_cov'][0],'b-*',label='$\\alpha = 3$')
plt.plot(ref_file['lamda_sc'],ref_file['rate_cov'][1],'g-*',label='$\\alpha = 3.5$')
plt.plot(ref_file['lamda_sc'],ref_file['rate_cov'][2],'r-*',label='$\\alpha = 4$')
plt.legend(loc='best')
plt.xlabel('$\lambda_{3}$')
plt.ylabel('rate coverage $\mathcal{R}_{\kappa_p}$')
plt.show()
#--------------------------------------fig7a----------------------------------------------------------------
'''
'''
#--------------------------------------fig8a----------------------------------------------------------------
ref_file = np.load('fig8a.npz')
plt.plot(ref_file['r_th2_ary']/10**6,ref_file['rate_cov'][0],'b-*',label='fair allocation $\kappa_{p}$')
plt.plot(ref_file['r_th2_ary']/10**6,ref_file['rate_cov'][1],'g-*',label='community ignorant')
plt.legend(loc='best')
plt.xlabel('$\\rho_{2} (\\times 10^{6})$')
plt.ylabel('rate coverage $\mathcal{R}$')
plt.show()
#--------------------------------------fig8a----------------------------------------------------------------
'''
'''
#--------------------------------------fig9a----------------------------------------------------------------
ref_file = np.load('fig9a.npz')
plt.plot(ref_file['r_th2_ary']/10**6,ref_file['rate_cov'][0],'b-*',label='fair allocation $\mathcal{R}^{1}$')
plt.plot(ref_file['r_th2_ary']/10**6,ref_file['rate_cov'][1],'g-*',label='fair allocation $\mathcal{R}^{2}$')
plt.plot(ref_file['r_th2_ary']/10**6,ref_file['rate_cov'][2],'r-*',label='comm ignorant $\mathcal{R}^{1}$')
plt.plot(ref_file['r_th2_ary']/10**6,ref_file['rate_cov'][3],'y-*',label='comm ignorant $\mathcal{R}^{2}$')
plt.legend(loc='best')
plt.xlabel('$\\rho_{2} (\\times 10^{6})$')
plt.ylabel('rate coverage')
plt.show()
#--------------------------------------fig9a----------------------------------------------------------------
'''

