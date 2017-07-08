from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
import matplotlib.pyplot as plt
import rate_coverage_v1 as r_c_v1

########################################################################################################################
           # define the system
########################################################################################################################

lamda_u1 = 0
lamda_u2 = 25
bw = 10*10**6
kappa = 0.5

#                     0       1            2              3        4          5          6
#                   [M-rat, K-tier, openaccess_flag, Pow(dBm), density,     alpha,    Bias(dB)]
m_1_k_1 =  np.array([1,       1,           True,          53,        1,         4,         0],float)
m_1_k_2 =  np.array([1,       2,           True,          33,        5,         4,         0],float)
m_1_k_3 = np.array([1,        3,           True,          23,       10,         4,         0],float)
mk_mat1 = np.array([ m_1_k_1, m_1_k_2, m_1_k_3])
lamda_u = 50
rate_threshold1 = ((np.array(range(0,22,1),float)/10)*10**6)
rate_threshold1_v1 = np.arange(0,2.3,0.1)*10**6

########################################################################################################################
           # rate coverage for two community of users
########################################################################################################################

def rate_cov_two_community(mk_matrix,lamda_u1,lamda_u2,th_u1,th_u2,bw,kappa):
    rate_cov1 = r_c_v1.rate_coverage_meanload(mk_matrix,lamda_u1,th_u1,kappa*bw)
    rate_cov2 = r_c_v1.rate_coverage_meanload(mk_matrix,lamda_u2,th_u2,(1-kappa)*bw)
    #print(rate_cov1,rate_cov2)
    rate_cov = 1 if (lamda_u1==lamda_u2==0) else (lamda_u1/(lamda_u1+lamda_u2))*rate_cov1 + (lamda_u2/(lamda_u1+lamda_u2))*rate_cov2
    return (rate_cov)
########################################################################################################################

'''
cr_r = rate_cov_two_community(mk_mat1,lamda_u1,lamda_u2,1*10**6,1*10**6,bw,kappa)
print(cr_r)
'''

'''
k = r_c_v1.rate_coverage_meanload(mk_mat1,lamda_u1,0*10**6,kappa*bw)
print(k)
'''

########################################################################################################################
           # rate coverage for two community of users - interference limited
########################################################################################################################

def rate_cov_two_community_int_lmt(mk_matrix,lamda_u1,lamda_u2,th_u1,th_u2,bw,kappa):
    rate_cov1 = r_c_v1.rate_coverage_meanload_int_lmt(mk_matrix,lamda_u1,th_u1,kappa*bw)
    rate_cov2 = r_c_v1.rate_coverage_meanload_int_lmt(mk_matrix,lamda_u2,th_u2,(1-kappa)*bw)
    #print(rate_cov1,rate_cov2)
    rate_cov = 1 if (lamda_u1==lamda_u2==0) else (lamda_u1/(lamda_u1+lamda_u2))*rate_cov1 + (lamda_u2/(lamda_u1+lamda_u2))*rate_cov2
    return (rate_cov)
########################################################################################################################

########################################################################################################################
           # rate coverage for two community of users - community ignorant
########################################################################################################################

def rate_cov_two_community_com_ignorant(mk_matrix,lamda_u1,lamda_u2,th_u1,th_u2,bw,kappa):
    rate_cov1 = r_c_v1.rate_coverage_meanload_int_lmt(mk_matrix,lamda_u1,th_u1,kappa*bw)
    rate_cov2 = r_c_v1.rate_coverage_meanload_int_lmt(mk_matrix,lamda_u2,th_u2,(1-kappa)*bw)
    #print(rate_cov1,rate_cov2)
    rate_cov = 1 if (lamda_u1==lamda_u2==0) else (lamda_u1/(lamda_u1+lamda_u2))*rate_cov1 + (lamda_u2/(lamda_u1+lamda_u2))*rate_cov2
    return (rate_cov)
########################################################################################################################
