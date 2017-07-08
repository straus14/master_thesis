from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
import matplotlib.pyplot as plt
import random_sleeping as rs
import pylab as pl
########################################################################################################################
 # system given is, \rho_{1} & \rho_{2} with threshold \e_{1} and \e_{2} - find optimum or minimum density of BS
########################################################################################################################
#-----------------------------------------------------------------------------------------------------------------------
# fairness optimization with common bandwidth approach - returns the denstiy of small cell required for the given rate coverage in the specific system setting
#-----------------------------------------------------------------------------------------------------------------------
def q_opt_rate_cov_fairness_v14(k_matrix,alpha,rate_th1,rate_th2,lamda_u1,lamda_u2,bw,tar_cov1,tar_cov2):
    k_mat1 = k_matrix.copy()
    rate_cov1 = 0
    rate_cov2 = 0
    kappa = 0
    q = 0
    while(rate_cov1 < tar_cov1 or rate_cov2 < tar_cov2 ):
        while(kappa <= 1):
            #print(type(q))
            rate_cov1 = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th1,lamda_u1,q, kappa*bw)
            rate_cov1 = round(rate_cov1,2)
            if ( rate_cov1 >= tar_cov1):
                break
            else:
                kappa = round(kappa + 0.1, 2)
        if ( kappa > 1):
            q = round(q+0.5,2)
            kappa = 0
        else:
            rate_cov2 = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th2,lamda_u2,q,(1-kappa)*bw)
            rate_cov2 = round(rate_cov2,2)
            if(rate_cov2 >= tar_cov2 and rate_cov1 >= tar_cov1):
                break
            else:
                q = round(q + 0.5,2)
                kappa = 0
    return (q)
#-----------------------------------------------------------------------------------------------------------------------

k1 = np.array([1,         43,             1,          10*10**6],float)
k2 = np.array([2,         38,             5,          10*10**6],float)
k3 = np.array([3,         21,             1,          10*10**6],float)
k_mat1 = np.array([k1, k2, k3],float)
alpha = 4

rate_th1 = 1 * 10**6
rate_th2 = 3 * 10**6
lamda_u1 = 25
lamda_u2 = 35
tar_cov1 = 0.5
tar_cov2 = 0.5
#tar_cov = 0.5
bw = 10*10**6
#-----------------------------------------------------------------------------------------------------------------------
'''
lamda_sc_fair = q_opt_rate_cov_fairness_v14(k_mat1,alpha,rate_th1,rate_th2,lamda_u1,lamda_u2,bw,tar_cov1,tar_cov2)
print(lamda_sc_fair)
'''

