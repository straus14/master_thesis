__author__ = 'straus14'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import pylab as pl

#-----------------------------------------------------------------------------------------------------------------------
          #figures for chapter 4 - UPDATED!!!
#-----------------------------------------------------------------------------------------------------------------------
fontP = FontProperties()
fontP.set_size('small')
'''
#-----------------------------------------------------------------------------------------------------------------------
          #variation of \mathcal_{R} vs B_{sc}^{(2)} for constant B_{sc}^{(1)} fig1
#-----------------------------------------------------------------------------------------------------------------------
ref_file = np.load('fig1.npz')
print(ref_file.files)
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][0],'b.-',label= '$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][0])+ 'dB' )
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][1],'gx-',label='$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][1])+ 'dB' )
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][2],'ys-',label='$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][2])+ 'dB')
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][3],'rh-',label='$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][3])+ 'dB')
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][4],'m*-',label='$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][4])+ 'dB')
plt.legend(loc=4,prop=fontP)
plt.xlabel('Association Bias for RAT-2 APs of $\mathcal{C}_{2}$ users - $B_{2}^{(2)}$  (dB)')
plt.ylabel('rate coverage ($\mathcal{R}$)')
plt.show()
#-----------------------------------------------------------------------------------------------------------------------
'''
'''
#-----------------------------------------------------------------------------------------------------------------------
          #fig2
#-----------------------------------------------------------------------------------------------------------------------
ref_file = np.load('fig2.npz')
print(ref_file.files)
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][0],'b.-',label= '$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][0])+ 'dB' )
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][1],'gx-',label='$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][1])+ 'dB' )
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][2],'ys-',label='$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][2])+ 'dB')
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][3],'rh-',label='$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][3])+ 'dB')
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][4],'m*-',label='$B_{2}^{(1)}$ = '+ str(ref_file['bias_c1'][4])+ 'dB')
plt.legend(loc=4,prop=fontP)
plt.xlabel('Association Bias for RAT-2 APs of $\mathcal{C}_{2}$ users - $B_{2}^{(2)}$  (dB)')
plt.ylabel('rate coverage ($\mathcal{R}$)')
plt.show()
#-----------------------------------------------------------------------------------------------------------------------
'''
'''
#-----------------------------------------------------------------------------------------------------------------------
          #fig3
#-----------------------------------------------------------------------------------------------------------------------
ref_file = np.load('fig3.npz')
print(ref_file.files)
ax = plt.subplot(111)
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][0],'b.-',label= '$\lambda_{2}$ = '+ str(ref_file['density_sc'][0])+'$\lambda_{1}$')
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][1],'gx-',label='$\lambda_{2}$ = '+ str(ref_file['density_sc'][1])+'$\lambda_{1}$')
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][2],'ys-',label='$\lambda_{2}$ = '+ str(ref_file['density_sc'][2])+'$\lambda_{1}$')
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][3],'rh-',label='$\lambda_{2}$ = '+ str(ref_file['density_sc'][3])+'$\lambda_{1}$')
plt.legend(loc='best',prop=fontP)
plt.xlabel('Association Bias for RAT-2 APs of $\mathcal{C}_{2}$ users - $B_{2}^{(2)}$  (dB)')
plt.ylabel('rate coverage ($\mathcal{R}$)')
plt.show()
#-----------------------------------------------------------------------------------------------------------------------
'''
'''
#-----------------------------------------------------------------------------------------------------------------------
          #fig4
#-----------------------------------------------------------------------------------------------------------------------
ref_file = np.load('fig4.npz')
print(ref_file.files)
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][0],'b.-',label='$\lambda_{u2}$ =' +str(ref_file['lamda_u2_ary'][0]))
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][1],'gx-',label='$\lambda_{u2}$ =' +str(ref_file['lamda_u2_ary'][1]))
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][2],'ys-',label='$\lambda_{u2}$ =' +str(ref_file['lamda_u2_ary'][2]))
plt.plot(ref_file['bias_c2'],ref_file['result_rate_cov'][3],'rh-',label='$\lambda_{u2}$ =' +str(ref_file['lamda_u2_ary'][3]))
plt.legend(loc=4,prop=fontP)
plt.xlabel('Association Bias for RAT-2 APs of $\mathcal{C}_{2}$ users - $B_{2}^{(2)}$  (dB)')
plt.ylabel('rate coverage ($\mathcal{R}$)')
plt.show()
#-----------------------------------------------------------------------------------------------------------------------
'''