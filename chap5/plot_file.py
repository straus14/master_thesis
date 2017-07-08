from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
import matplotlib.pyplot as plt
import pylab as pl
import two_RAT_rate_cov as offload1
########################################################################################################################
          # rate coverage for two RAT system two community system - figures
########################################################################################################################

########################################################################################################################
   # system model
########################################################################################################################
#                 0            1             2               3           4                 5
#                [K-tier,    Pow(dBm),     density,         BW,         Bias_C1(dB),   Bias_C2(dB) ]
k1 =    np.array([1,         53,             1,          10*10**6,        0,              0],float)
k2 =    np.array([2,         23,             10,         10*10**6,        0,             0],float)
k_mat1 = np.array([k1, k2],float)

'''
########################################################################################################################
                                    # chap4 - fig1
########################################################################################################################

alpha = 4
lamda_u1 = 25
lamda_u2 = 75
rate_th1 = 0.5*10**6
rate_th2 = 2*10**6

bias_c1 = pl.frange(10,30,5)
bias_c2 = pl.frange(0,80,3)
result_rate_cov = np.zeros((len(bias_c1),len(bias_c2)))
for idx1,val1 in enumerate(bias_c1):
    k_mat1[1][4] = val1
    for idx2,val2 in enumerate(bias_c2):
        k_mat1[1][5]=val2
        result_rate_cov[idx1][idx2] = offload1.rate_cov_two_rat_two_com(k_mat1,alpha,lamda_u1,lamda_u2,rate_th1,rate_th2)
np.savez('/home/straus14/hiwi_python/task_1/bias_community_two_RAT/fig1',bias_c1=bias_c1,bias_c2=bias_c2,result_rate_cov=result_rate_cov)
########################################################################################################################
'''
'''
########################################################################################################################
                                    # chap4 - fig2
########################################################################################################################

alpha = 4
lamda_u1 = 25
lamda_u2 = 25
rate_th1 = 0.5*10**6
rate_th2 = 4*10**6

bias_c1 = pl.frange(0,20,5)
bias_c2 = pl.frange(0,80,3)
result_rate_cov = np.zeros((len(bias_c1),len(bias_c2)))
for idx1,val1 in enumerate(bias_c1):
    k_mat1[1][4] = val1
    for idx2,val2 in enumerate(bias_c2):
        k_mat1[1][5]=val2
        result_rate_cov[idx1][idx2] = offload1.rate_cov_two_rat_two_com(k_mat1,alpha,lamda_u1,lamda_u2,rate_th1,rate_th2)
np.savez('/home/straus14/hiwi_python/task_1/bias_community_two_RAT/fig2',bias_c1=bias_c1,bias_c2=bias_c2,result_rate_cov=result_rate_cov)
########################################################################################################################
'''
'''
########################################################################################################################
                                    # chap4 - fig3
########################################################################################################################

alpha = 4
lamda_u1 = 25
lamda_u2 = 75
rate_th1 = 0.5*10**6
rate_th2 = 2*10**6

density_sc = pl.frange(5,20,5)
bias_c2 = pl.frange(0,40,3)
result_rate_cov = np.zeros((len(density_sc),len(bias_c2)))
for idx1,val1 in enumerate(density_sc):
    k_mat1[1][2]=val1
    for idx2,val2 in enumerate(bias_c2):
        k_mat1[1][5] = val2
        result_rate_cov[idx1][idx2] = offload1.rate_cov_two_rat_two_com(k_mat1,alpha,lamda_u1,lamda_u2,rate_th1,rate_th2)
np.savez('/home/straus14/hiwi_python/task_1/bias_community_two_RAT/fig3',density_sc=density_sc,bias_c2=bias_c2,result_rate_cov=result_rate_cov)

########################################################################################################################
'''
'''
########################################################################################################################
                                    # chap4 - fig4
########################################################################################################################

alpha = 4
lamda_u1 = 25
lamda_u2 = 75
rate_th1 = 0.5*10**6
rate_th2 = 2*10**6

lamda_u2_ary = pl.frange(25,100,25)
bias_c2 = pl.frange(0,40,3)
result_rate_cov = np.zeros((len(lamda_u2_ary),len(bias_c2)))
for idx1,val1 in enumerate(lamda_u2_ary):
    lamda_u2=val1
    for idx2,val2 in enumerate(bias_c2):
        k_mat1[1][5] = val2
        result_rate_cov[idx1][idx2] = offload1.rate_cov_two_rat_two_com(k_mat1,alpha,lamda_u1,lamda_u2,rate_th1,rate_th2)
np.savez('/home/straus14/hiwi_python/task_1/bias_community_two_RAT/fig4',lamda_u2_ary=lamda_u2_ary,bias_c2=bias_c2,result_rate_cov=result_rate_cov)

########################################################################################################################
'''