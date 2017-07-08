from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
import matplotlib.pyplot as plt
import rate_coverage_v1 as r_c_v1
import pylab as pl
import two_community_v1 as r_c_v2
########################################################################################################################


########################################################################################################################
              # define the system
########################################################################################################################

#                     0       1            2              3        4          5          6
#                   [M-rat, K-tier, openaccess_flag, Pow(dBm), density,     alpha,    Bias(dB)]
m_1_k_1 =  np.array([1,       1,           True,          53,        1,         4,         0],float)
m_1_k_2 =  np.array([1,       2,           True,          33,        5,         4,         0],float)
m_1_k_3 = np.array([1,        3,           True,          23,       10,         4,         0],float)
mk_mat1 = np.array([ m_1_k_1, m_1_k_2, m_1_k_3])

'''
########################################################################################################################
           # chap3-fig1  - rate coverage vs rate threshold two community - for fixed user density
########################################################################################################################
lamda_u1 = 38
lamda_u2 = 12
bw = 10*10**6
rate_th1 = pl.frange(0,6,0.5)*10**6
rate_th2 = pl.frange(0,6,0.5)*10**6
rate_cov_result_fig1 = np.zeros((len(rate_th1),len(rate_th2)))

for idx1,val1 in enumerate(rate_th1):
    r_th1 = val1
    for idx2,val2 in enumerate(rate_th2):
        r_th2 = val2
        kappa = 0 if (r_th1==r_th2==0 or lamda_u1==lamda_u2==0) else (lamda_u1*r_th1)/(lamda_u1*r_th1 + lamda_u2 * r_th2)
        rate_cov_result_fig1[idx1][idx2] = r_c_v2.rate_cov_two_community(mk_mat1,lamda_u1,lamda_u2,r_th1,r_th2,bw,kappa)

np.savez('/home/straus14/hiwi_python/task_1/Community_of_Users/fig1b',rate_th1=rate_th1,rate_th2=rate_th2,rate_cov_result_fig1b=rate_cov_result_fig1)


#########################################################################################################################
'''
'''
########################################################################################################################
           # chap3-fig2  - rate coverage vs rate threshold two community - for fixed user density
########################################################################################################################
bw = 10*10**6
rate_th1 = 0.5*10**6
rate_th2 = 4*10**6
lamda_u1_ary = pl.frange(0,50,5)
lamda_u2_ary = pl.frange(0,50,5)
rate_cov_result_fig2 = np.zeros((len(lamda_u1_ary),len(lamda_u1_ary)))

for idx1,val1 in enumerate(lamda_u1_ary):
    lamda_u1 = val1
    for idx2,val2 in enumerate(lamda_u2_ary):
        lamda_u2 = val2
        kappa = 0 if (rate_th1==rate_th2==0 or lamda_u1==lamda_u2==0) else (lamda_u1*rate_th1)/(lamda_u1*rate_th1 + lamda_u2 * rate_th2)
        rate_cov_result_fig2[idx1][idx2] = r_c_v2.rate_cov_two_community(mk_mat1,lamda_u1,lamda_u2,rate_th1,rate_th2,bw,kappa)

np.savez('/home/straus14/hiwi_python/task_1/Community_of_Users/fig2b',lamda_u1_ary=lamda_u1_ary,lamda_u2_ary=lamda_u2_ary,rate_cov_result_fig2b=rate_cov_result_fig2)


#########################################################################################################################
'''
'''
########################################################################################################################
           # chap3-fig3  - rate coverage vs rate threshold two community - maximum of rate
########################################################################################################################
lamda_u1 = 38
lamda_u2 = 12
bw = 10*10**6
rate_th1 = pl.frange(0,6,0.5)*10**6
rate_th2 = pl.frange(0,6,0.5)*10**6
rate_cov_result_fig3 = np.zeros((len(rate_th1),len(rate_th2)))

for idx1,val1 in enumerate(rate_th1):
    r_th1 = val1
    for idx2,val2 in enumerate(rate_th2):
        r_th2 = val2
        r_th = max(val1,val2)
        kappa = 0 if (r_th1==r_th2==0 or lamda_u1==lamda_u2==0) else (lamda_u1*r_th1)/(lamda_u1*r_th1 + lamda_u2 * r_th2)
        rate_cov_result_fig3[idx1][idx2] = r_c_v2.rate_cov_two_community(mk_mat1,lamda_u1,lamda_u2,r_th,r_th,bw,kappa)

np.savez('/home/straus14/hiwi_python/task_1/Community_of_Users/fig3b',rate_th1=rate_th1,rate_th2=rate_th2,rate_cov_result_fig3b=rate_cov_result_fig3)
#########################################################################################################################
'''
'''
########################################################################################################################
           # chap3-fig4  - variation of rate coverage with community
########################################################################################################################
kappa_array = pl.frange(0,1,0.05)
lamda_u1 = 300
lamda_u2 = 200
r_th1 = 2*10**6
r_th2 = 4*10**6
bw = 50 * 10**6
rate_cov1 = np.zeros(len(kappa_array))
rate_cov2 = np.zeros(len(kappa_array))
rate_cov = np.zeros(len(kappa_array))
for idx,val in enumerate(kappa_array):
    rate_cov1[idx] = 0 if val==0 else r_c_v1.rate_coverage_meanload(mk_mat1,lamda_u1,r_th1,val*bw)
    rate_cov2[idx] = 0 if (1-val)==0 else r_c_v1.rate_coverage_meanload(mk_mat1,lamda_u2,r_th2,(1-val)*bw)
    rate_cov[idx] = (lamda_u1/(lamda_u1+lamda_u2))*rate_cov1[idx] + (lamda_u2/(lamda_u1+lamda_u2))*rate_cov2[idx]
np.savez('/home/straus14/hiwi_python/task_1/Community_of_Users/fig5a',kappa_array=kappa_array,rate_cov1=rate_cov1,rate_cov2=rate_cov2,rate_cov=rate_cov)
########################################################################################################################
'''
'''
########################################################################################################################
           # chap3-fig6  - rate coverage variation with base station density
########################################################################################################################
lamda_u1 = 25
lamda_u2 = 25
r_th1 = 0.5*10**6
r_th2 = 4*10**6
bw = 10 * 10**6
kappa = (lamda_u1*r_th1)/(lamda_u1*r_th1 + lamda_u2 * r_th2)
lamda_sc = range(0,200,10)
lamda_mbs = range(0,15,5)
rate_cov = np.zeros((len(lamda_mbs),len(lamda_sc)))
for idx1,val1 in enumerate(lamda_mbs):
    mk_mat1[0][4] = val1
    for idx2,val2 in enumerate(lamda_sc):
        mk_mat1[2][4] = val2
        rate_cov[idx1][idx2] = r_c_v2.rate_cov_two_community(mk_mat1,lamda_u1,lamda_u2,r_th1,r_th2,bw,kappa)

np.savez('/home/straus14/hiwi_python/task_1/Community_of_Users/fig6a',lamda_mbs=lamda_mbs,lamda_sc=lamda_sc,rate_cov=rate_cov)
########################################################################################################################
'''
'''
########################################################################################################################
           # chap3-fig7  - rate coverage variation with base station density with various alphas
########################################################################################################################
lamda_u1 = 25
lamda_u2 = 25
r_th1 = 0.5*10**6
r_th2 = 4*10**6
bw = 10 * 10**6
kappa = (lamda_u1*r_th1)/(lamda_u1*r_th1 + lamda_u2 * r_th2)
alpha_ary = pl.frange(3,4,0.5)
lamda_sc = range(0,105,10)
rate_cov = np.zeros((len(alpha_ary),len(lamda_sc)))
for idx1,val1 in enumerate(alpha_ary):
    mk_mat1[0][5] = val1
    mk_mat1[1][5] = val1
    mk_mat1[2][5] = val1
    for idx2,val2 in enumerate(lamda_sc):
        mk_mat1[2][4] = val2
        rate_cov[idx1][idx2] = r_c_v2.rate_cov_two_community(mk_mat1,lamda_u1,lamda_u2,r_th1,r_th2,bw,kappa)

np.savez('/home/straus14/hiwi_python/task_1/Community_of_Users/fig7a',alpha_ary=alpha_ary,lamda_sc=lamda_sc,rate_cov=rate_cov)
########################################################################################################################
'''
'''
########################################################################################################################
           # chap3-fig8  - fair allocation, community ignorant, proposition-2
########################################################################################################################
lamda_u1 = 30
lamda_u2 = 15
r_th1 = 0.5*10**6
r_th2 = 4*10**6
bw = 10 * 10**6
r_th2_ary = pl.frange(0,6,0.5)*10**6
rate_cov = np.zeros((2,len(r_th2_ary)))
for idx,val in enumerate(r_th2_ary):
    r_th2 = val
    kappa_fair = (lamda_u1*r_th1)/(lamda_u1*r_th1 + lamda_u2 * r_th2)
    kappa_prop2 = (lamda_u1)/(lamda_u1 + lamda_u2)
    rate_cov[0][idx] = r_c_v2.rate_cov_two_community(mk_mat1,lamda_u1,lamda_u2,r_th1,r_th2,bw,kappa_fair)
    rate_cov[1][idx] = r_c_v2.rate_cov_two_community(mk_mat1,lamda_u1,lamda_u2,r_th1,r_th2,bw,kappa_prop2)

np.savez('/home/straus14/hiwi_python/task_1/Community_of_Users/fig8a',r_th2_ary=r_th2_ary,rate_cov=rate_cov)
########################################################################################################################
'''
'''
########################################################################################################################
           # chap3-fig9  - fair allocation, community ignorant
########################################################################################################################
lamda_u1 = 30
lamda_u2 = 15
r_th1 = 0.5*10**6
r_th2 = 4*10**6
bw = 10 * 10**6
r_th2_ary = pl.frange(0,6,0.5)*10**6
rate_cov = np.zeros((4,len(r_th2_ary)))
for idx,val in enumerate(r_th2_ary):
    r_th2 = val
    kappa_fair = (lamda_u1*r_th1)/(lamda_u1*r_th1 + lamda_u2 * r_th2)
    kappa_ignorant = (lamda_u1)/(lamda_u1 + lamda_u2)
    rate_cov[0][idx] = r_c_v1.rate_coverage_meanload(mk_mat1,lamda_u1,r_th1,bw*kappa_fair)
    rate_cov[1][idx] = r_c_v1.rate_coverage_meanload(mk_mat1,lamda_u2,r_th2,bw*(1-kappa_fair))
    rate_cov[2][idx] = r_c_v1.rate_coverage_meanload(mk_mat1,lamda_u1,r_th1,bw*kappa_ignorant)
    rate_cov[3][idx] = r_c_v1.rate_coverage_meanload(mk_mat1,lamda_u2,r_th2,bw*(1-kappa_ignorant))

np.savez('/home/straus14/hiwi_python/task_1/Community_of_Users/fig9a',r_th2_ary=r_th2_ary,rate_cov=rate_cov)
########################################################################################################################
'''
