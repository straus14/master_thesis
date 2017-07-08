from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
import matplotlib.pyplot as plt
import pylab as pl
########################################################################################################################

'''
########################################################################################################################
#chapter 32 - fig2
########################################################################################################################
ref_file = np.load('fig2.npz')
plt.plot(ref_file['lamda_sc'],ref_file['EE_mat'][0]*10**3,'b*-',label='$\\rho_{_2}=1.5\\times10^{6}$')
plt.plot(ref_file['lamda_sc'],ref_file['EE_mat'][1]*10**3,'g--',label='$\\rho_{_2}=2\\times10^{6}$')
plt.plot(ref_file['lamda_sc'],ref_file['EE_mat'][2]*10**3,'r+-',label='$\\rho_{_2}=2.5\\times10^{6}$')
plt.plot(ref_file['lamda_sc'],ref_file['EE_mat'][3]*10**3,'yo-',label='$\\rho_{_2}=3\\times10^{6}$')
plt.legend(loc='best')
plt.xlabel('small cell BS density ($\lambda_{3}$)')
plt.ylabel('Energy Efficiency (Bits/s/Hz/W) $\\times 10^{-3}$')
plt.show()
########################################################################################################################
'''
'''
########################################################################################################################
#chapter 32 - fig3
########################################################################################################################
ref_file = np.load('fig3.npz')
k = ref_file['rate_th2_aray']
plt.plot(ref_file['rate_th2_aray']/10**6,ref_file['result_lamda_sc'][0],'gs-',label='$\epsilon_{2}$=0.5')
plt.plot(ref_file['rate_th2_aray']/10**6,ref_file['result_lamda_sc'][1],'b*-',label='$\epsilon_{2}$=0.4')
plt.plot(ref_file['rate_th2_aray']/10**6,ref_file['result_lamda_sc'][2],'ro-',label='$\epsilon_{2}$=0.3')
plt.legend(loc='best')
plt.xlabel('rate threshold for community-$\mathcal{C}_{2}$ ($\\rho_{_2}$) ($\\times 10^{6}$)')
plt.ylabel('required minimum small cell density $\lambda_{3}$')
plt.show()
print(ref_file.files)
########################################################################################################################
'''
'''
########################################################################################################################
#chapter 32 - fig4
########################################################################################################################
ref_file = np.load('fig4.npz')
plt.plot(ref_file['tar_cov2_ary'],ref_file['result_lamda_sc'][0],'gs-',label='$\\rho_{2}=2\\times10^{6}$')
plt.plot(ref_file['tar_cov2_ary'],ref_file['result_lamda_sc'][1],'b*-',label='$\\rho_{2}=3\\times10^{6}$')
plt.plot(ref_file['tar_cov2_ary'],ref_file['result_lamda_sc'][2],'ro-',label='$\\rho_{2}=4\\times10^{6}$')
plt.legend(loc='best')
plt.xlabel('required rate coverage for community $\mathcal{C}_{2}$($\epsilon_{2}$)')
plt.ylabel('required minimum small cell density $\lambda_{3}$')
plt.show()
print(ref_file.files)
print(ref_file['tar_cov2_ary'])
########################################################################################################################
'''
'''
########################################################################################################################
#chapter 32 - fig5
########################################################################################################################
ref_file = np.load('fig5.npz')
plt.plot(ref_file['lamda_u2_ary'],ref_file['result_lamda_sc'][0],'gs-',label='$\\rho_{2}=2\\times10^{6}$')
plt.plot(ref_file['lamda_u2_ary'],ref_file['result_lamda_sc'][1],'b*-',label='$\\rho_{2}=3\\times10^{6}$')
plt.plot(ref_file['lamda_u2_ary'],ref_file['result_lamda_sc'][2],'ro-',label='$\\rho_{2}=4\\times10^{6}$')
plt.legend(loc='best')
plt.xlabel('user density for community $\mathcal{C}_{2}$($\lambda_{u2}$)')
plt.ylabel('required minimum small cell density $\lambda_{3}$')
plt.show()
print(ref_file.files)
########################################################################################################################
'''
'''
########################################################################################################################
#chapter 32 - fig6
########################################################################################################################
ref_file = np.load('fig6.npz')
plt.plot(ref_file['lamda_u2_ary'],ref_file['result_lamda_sc'][0],'gs-',label='$\\epsilon_{2}=0.3$')
plt.plot(ref_file['lamda_u2_ary'],ref_file['result_lamda_sc'][1],'b*-',label='$\\epsilon_{2}=0.5$')
plt.plot(ref_file['lamda_u2_ary'],ref_file['result_lamda_sc'][2],'ro-',label='$\\epsilon_{2}=0.7$')
plt.legend(loc='best')
plt.xlabel('user density for community $\mathcal{C}_{2}$($\lambda_{u2}$)')
plt.ylabel('required minimum small cell density $\lambda_{3}$')
plt.show()
print(ref_file.files)
########################################################################################################################
'''
'''
########################################################################################################################
#chapter 32 - fig7
########################################################################################################################
ref_file = np.load('fig7.npz')
print(ref_file.files)
plt.plot(ref_file['q_on'],ref_file['ee_rs']*10**3,'r+-',label='random sleeping')
plt.plot(ref_file['q_on'],ref_file['ee_ss']*10**3,'b*-',label='strategic sleeping')
plt.legend(loc='best')
plt.ylabel('Energy Efficiency (Bits/s/Hz/W) ($\\times 10^{-3}$)')
plt.xlabel('fraction of SCBS awake (q_on)')
plt.show()
########################################################################################################################
'''
'''
########################################################################################################################
#chapter 32 - fig8
########################################################################################################################
ref_file = np.load('fig8.npz')
print(ref_file.files)
plt.plot(ref_file['rate_th_ary']/10**6,ref_file['rate_cov_ss_binary'],'b*-',label='binary')
plt.plot(ref_file['rate_th_ary']/10**6,ref_file['rate_cov_ss_uniform'],'r+-',label='uniform')
plt.legend(loc='best')
plt.xlabel('rate threshold ($\\times 10^{6}$)')
plt.ylabel('rate coverage')
plt.show()
########################################################################################################################
'''
'''
########################################################################################################################
#chapter 32 - fig9
########################################################################################################################
multiple_user_SC_density=np.array([  0.  ,   5.25,  23.75,  38.5 ,  51.75,  64.5 ,  77.5])
single_user_SC_density=np.array([ 0.  ,  13.5 ,  35.  ,  52.25,  68.  ,  83.25,  98.25])
rate_th_array = np.array([       0.0,   500000.0,  1000000.0,  1500000.0,  2000000.0,  2500000.0,  3000000.0,])/(10**6)
#plt.figure(figsize=(3.5,2.5))
#matplotlib.rcParams.update({'font.size': 7})
plt.plot(rate_th_array,multiple_user_SC_density,'bo-',rate_th_array,single_user_SC_density,'rs-')
#plt.xlabel('Two community:$\\rho^{1}$ ($\\rho^{2}$ = 2$\\times\\rho^{1}$), Single Community: $max(\\rho^{1},\\rho^{2})$')
plt.xlabel('$\\rho^{1} (\\times 10^{6})$')
plt.ylabel('small cell density $\lambda_{3}$')
#plt.legend(['two community system','single community system'],loc=4)
plt.legend(['two community system','single community system'],'upper left')
#plt.tight_layout()
#plt.savefig('example_2.pgf')
#plt.legend(loc='best')
plt.show()
########################################################################################################################
'''