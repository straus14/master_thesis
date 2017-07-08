from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
import matplotlib.pyplot as plt
import pylab as pl
########################################################################################################################
          # rate coverage for two RAT system two community system
########################################################################################################################

def rate_cov_two_rat_two_com(mk_matrix,alpha,lamda_u1,lamda_u2,rate_th1,rate_th2):
    mk_mat = mk_matrix.copy()
    mk_mat_size = mk_mat.shape
    mk_mat[:,1] = gb.dbm2pow(mk_mat[:,1])
    mk_mat[:,4] = gb.db2pow(mk_mat[:,4])       # bias for C1 users
    mk_mat[:,5] = gb.db2pow(mk_mat[:,5])       # bias for C2 users
    density = mk_mat[:,2].copy()
    power = mk_mat[:,1].copy()
    bias_c1 = mk_mat[:,4].copy()
    bias_c2 = mk_mat[:,5].copy()
    # initialize the needed values
    area_c1 = np.zeros(mk_mat_size[0])         # association probability for community1 users
    area_c2 = np.zeros(mk_mat_size[0])         # association probability for community2 users
    N_ij_c1 = np.zeros(mk_mat_size[0])         # number of users in for C1 users
    N_ij_c2 = np.zeros(mk_mat_size[0])         # number of users in for C2 users
    threshold_c1 = np.zeros(mk_mat_size[0])    # threshold for C1 users
    threshold_c2 = np.zeros(mk_mat_size[0])    # threshold for C2 users
    bw_c1 = np.zeros(mk_mat_size[0])           # bandwidth for C1 users
    bw_c2 = np.zeros(mk_mat_size[0])           # bandwidth for C2 users
    kappa = np.zeros(mk_mat_size[0])           # fraction of BW for C1 users
    rate_cov_c1 = np.zeros(mk_mat_size[0])     # rate coverage for C1 users
    rate_cov_c2 = np.zeros(mk_mat_size[0])     # rate coverage for C2 users
    temp_sum_c1 = mk_mat[:,2] * (mk_mat[:,1]* mk_mat[:,4])**(2/alpha)        #temporary sum for denominator of G function of C1
    temp_sum_c2 = mk_mat[:,2] * (mk_mat[:,1]* mk_mat[:,5])**(2/alpha)        #temporary sum for denominator of G function of C2
    #fill the values in arrays
    for i in range(mk_mat_size[0]):
        area_c1[i] = temp_sum_c1[i]/(sum(temp_sum_c1))                # area of C1
        area_c2[i] = temp_sum_c2[i]/(sum(temp_sum_c2))                # area of C2
        #kappa[i] = (lamda_u1*rate_th1*area_c1[i])/((lamda_u1*rate_th1*area_c1[i])+(lamda_u2*rate_th2*area_c2[i]))
    N_ij_c1 = 1 + 1.28*lamda_u1*(area_c1/mk_mat[:,2])
    N_ij_c2 = 1 + 1.28*lamda_u2*(area_c2/mk_mat[:,2])
    #N_ij_c1 =  1 + lamda_u1*(area_c1/mk_mat[:,2])
    #N_ij_c2 = 1 + lamda_u2*(area_c2/mk_mat[:,2])
    kappa = 0 if (rate_th1==0 and rate_th2==0) else (lamda_u1*rate_th1*area_c1)/((lamda_u1*rate_th1*area_c1)+(lamda_u2*rate_th2*area_c2))
    #kappa = np.array([0,0])
    bw_c1 = kappa * mk_mat[:,3]
    bw_c2 = (1-kappa) * mk_mat[:,3]
    for i in range(mk_mat_size[0]):
        threshold_c1[i] = 0 if (bw_c1[i]==0) else 2**((rate_th1*N_ij_c1[i])/(bw_c1[i])) - 1
        threshold_c2[i] = 0 if (bw_c2[i]==0) else 2**((rate_th2*N_ij_c2[i])/(bw_c2[i])) - 1
    #threshold_c1 = 2**((rate_th1*N_ij_c1)/(bw_c1)) - 1
    #threshold_c2 = 2**((rate_th2*N_ij_c2)/(bw_c2)) - 1
    # finding rate coverage for community 1 users
    for i in range(mk_mat_size[0]):
        z_func_c1 = 0 if (threshold_c1[i]==0) else (threshold_c1[i])**(2/alpha) * mp.quad(lambda u: (1/(1+u**(alpha/2))), [(1/threshold_c1[i])**(2/alpha), mp.inf])
        sec_term_denom_c1 = sum(temp_sum_c1)/((power[i]*bias_c1[i])**(2/alpha))
        rate_cov_c1[i] = density[i]/(density[i]*z_func_c1 + sec_term_denom_c1)
        z_func_c2 = 0 if (threshold_c2[i]==0) else (threshold_c2[i])**(2/alpha) * mp.quad(lambda u: (1/(1+u**(alpha/2))), [(1/threshold_c2[i])**(2/alpha), mp.inf])
        sec_term_denom_c2 = sum(temp_sum_c2)/((power[i]*bias_c2[i])**(2/alpha))
        rate_cov_c2[i] = density[i]/(density[i]*z_func_c2 + sec_term_denom_c2)
    total_rate_cov = (lamda_u1/(lamda_u1+lamda_u2)) * (sum(rate_cov_c1)) + (lamda_u2/(lamda_u1+lamda_u2)) * (sum(rate_cov_c2))
    ase = (lamda_u1*rate_th1*sum(rate_cov_c1) + lamda_u2*rate_th2*sum(rate_cov_c2))/(10*10**6)
    #print(threshold_c1)
    #print(threshold_c2)
    #print(rate_cov_c1,sum(rate_cov_c1))
    #print(rate_cov_c2,sum(rate_cov_c2))
    #print(kappa)
    #print(rate_cov_c1,sum(rate_cov_c1))
    #print(rate_cov_c2,sum(rate_cov_c2))
    #print(ase)
    #print(N_ij_c1)
    #print(N_ij_c2)
    return (total_rate_cov)

























'''
########################################################################################################################
   # system model
########################################################################################################################
#                 0            1             2               3           4                 5
#                [K-tier,    Pow(dBm),     density,         BW,         Bias_C1(dB),   Bias_C2(dB) ]
k1 =    np.array([1,         43,             1,          10*10**6,        0,              0],float)
k2 =    np.array([2,         21,             5,          10*10**6,       15,              20],float)
k_mat1 = np.array([k1, k2],float)

alpha = 4
lamda_u1 = 50
lamda_u2 = 50
rate_th1 = 0.5*10**6
rate_th2 = 8*10**6

k = rate_cov_two_rat_two_com(k_mat1,alpha,lamda_u1,lamda_u2,rate_th1,rate_th2)
print(k)
'''