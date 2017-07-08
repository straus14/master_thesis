#random switching
from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
from scipy.integrate import dblquad,quad
import matplotlib.pyplot as plt
#######################################################################################################################
#######################################################################################################################

#area calculation for tiers with equal pathloss exponent
def A_k(density,power,alpha,idx):
    area_k = density[idx] / (sum(density*(power/power[idx])**(2/alpha)))
    return area_k

#-----------------------------------------------------------------------------------------------------------------------

def rate_cov_strate_sleeping_v2(k_matrix,alpha,rate_th1,lamda_u1,bw,prob_actv):
    # define the distribution of the activity
    first_int = prob_actv                                     # correspond to the integral in first term
    second_int = 0                                            # correspond to the integral in second term
    #---------------------------- calibrated -----------------------------------------
    expected_activity = prob_actv                             # expected value of a      # this changes for optimization hance should be calibratable
    expected_strategic_function = prob_actv                           # expected value of s
    #---------------------------- calibrated -----------------------------------------
    noise_var = 1
    # preprocessing - define empty elements and get information about the inputs
    k_mat = k_matrix.copy()
    num_tiers = k_mat.shape[0]
    density_org = k_mat[:,2].copy()
    power = gb.db2pow(k_mat[:,1])
    density_update = np.array(density_org*([1]*(num_tiers-1)+[expected_strategic_function]))
    # define necessary values
    area_org = np.zeros(num_tiers,float)                      # original association probability
    area_sc_update = np.zeros(num_tiers,float)                # association probability of disconnected cell
    N_k_u1 = np.zeros(num_tiers,float)                        #number of users in tier K BS
    N_k_sc = np.zeros(num_tiers,float)                        #number of users in tier K BS
    N_k_total = np.zeros(num_tiers,float)
    threshold_u1 = np.zeros(num_tiers,float)                  #threshold for users in tier K BS
    t_func_main = np.zeros(num_tiers,float)
    t_func_sc = np.zeros(num_tiers,float)
    for i in range(num_tiers):
        area_org[i] = A_k(density_org,power,alpha,i)
        area_sc_update[i] = A_k(density_update,power,alpha,i)
    N_k_u1 = 1 + 1.28*lamda_u1*(area_org/density_org)
    #N_k_sc = expected_activity*(1-expected_strategic_function)*density_org[-1]*N_k_u1[-1]*(area_sc_update/density_update)   #for binary optimization
    N_k_sc = 0                              # set N_k_sc = 0 for binary activity modeling
    N_k_total = N_k_u1 + N_k_sc
    threshold_u1 = 2**((rate_th1/bw)*N_k_total) -1
    for i in range(num_tiers):
        first_exp_term = -(threshold_u1[i]*noise_var/power[i])
        z_term = 0 if (threshold_u1[i]==0) else (threshold_u1[i])**(2/alpha) *mp.quad(lambda u: 1/ (1+u**(alpha/2)),[(1/threshold_u1[i])**(2/alpha),mp.inf])
        second_exp_term = -mp.pi*z_term* sum(density_update*(power/power[i])**(2/alpha))
        third_exp_term_main = -mp.pi * sum(density_org*(power/power[i])**(2/alpha))
        third_exp_term_sc = -mp.pi * sum(density_update*(power/power[i])**(2/alpha))
        t_func_main[i] = mp.quad(lambda y: y * mp.exp(first_exp_term*y**alpha) * mp.exp(second_exp_term*y**2) * mp.exp(third_exp_term_main*y**2),[0,mp.inf])
        t_func_sc[i] = mp.quad(lambda y: y* mp.exp(first_exp_term*y**alpha) * mp.exp(second_exp_term*y**2) * mp.exp(third_exp_term_sc*y**2), [0,mp.inf])
    temp_second_sum = sum(2*mp.pi*density_update*t_func_sc)
    temp_third_sum = sum(2*mp.pi*density_org[0:-1]*t_func_main[0:-1])
    rate_coverage = (2*mp.pi*density_org[-1]/expected_activity)*(t_func_main[-1]*first_int + temp_second_sum*second_int) + temp_third_sum
    return (rate_coverage)










#-----------------------------------------------------------------------------------------------------------------------
#                 0            1             2               3
#                [K-tier,    Pow(dBm),     density,         BW ]
k1 =    np.array([1,         43,             1,          10*10**6],float)
k2 =    np.array([2,         38,             5,          10*10**6],float)
k3 =    np.array([3,         21,             10,         10*10**6],float)
k_mat1 = np.array([k1, k2, k3],float)
alpha = 4
rate_th1 = 0.5 * 10**6
rate_th2 = 1 * 10**6
lamda_u1 = 25
lamda_u2 = 25
lamda_user = 50
bw = 10*10**6
sc_act_prob = 0.1
'''
result1 = rate_cov_strate_sleeping_v2(k_mat1,alpha,rate_th1,lamda_u1,bw,sc_act_prob)
print(result1)'''
#-----------------------------------------------------------------------------------------------------------------------

def rate_cov_strate_sleeping_v3_bin(k_matrix,alpha,rate_th1,lamda_u1,bw,prob_actv):
    # define the distribution of the activity
    first_int = prob_actv                                     # correspond to the integral in first term
    second_int = 0                                            # correspond to the integral in second term
    #---------------------------- calibrated -----------------------------------------
    expected_activity = prob_actv                             # expected value of a      # this changes for optimization hance should be calibratable
    expected_strategic_function = prob_actv                           # expected value of s
    #---------------------------- calibrated -----------------------------------------
    noise_var = 1
    # preprocessing - define empty elements and get information about the inputs
    k_mat = k_matrix.copy()
    num_tiers = k_mat.shape[0]
    density_org = k_mat[:,2].copy()
    power = gb.db2pow(k_mat[:,1])
    density_update = np.array(density_org*([1]*(num_tiers-1)+[expected_strategic_function]))
    # define necessary values
    area_org = np.zeros(num_tiers,float)                      # original association probability
    area_sc_update = np.zeros(num_tiers,float)                # association probability of disconnected cell
    N_k_u1 = np.zeros(num_tiers,float)                        #number of users in tier K BS
    N_k_sc = np.zeros(num_tiers,float)                        #number of users in tier K BS
    N_k_total = np.zeros(num_tiers,float)
    threshold_u1 = np.zeros(num_tiers,float)                  #threshold for users in tier K BS
    t_func_main = np.zeros(num_tiers,float)
    t_func_sc = np.zeros(num_tiers,float)
    for i in range(num_tiers):
        area_org[i] = A_k(density_org,power,alpha,i)
        area_sc_update[i] = A_k(density_update,power,alpha,i)
    N_k_u1 = 1 + 1.28*lamda_u1*(area_org/density_org)
    N_k_sc = 0
    N_k_total = N_k_u1 + N_k_sc
    threshold_u1 = 2**((rate_th1/bw)*N_k_total) -1
    for i in range(num_tiers):
        first_exp_term = -(threshold_u1[i]*noise_var/power[i])
        z_term = 0 if (threshold_u1[i]==0) else (threshold_u1[i])**(2/alpha) *mp.quad(lambda u: 1/ (1+u**(alpha/2)),[(1/threshold_u1[i])**(2/alpha),mp.inf])
        second_exp_term = -mp.pi*z_term* sum(density_update*(power/power[i])**(2/alpha))
        third_exp_term_main = -mp.pi * sum(density_org*(power/power[i])**(2/alpha))
        third_exp_term_sc = -mp.pi * sum(density_update*(power/power[i])**(2/alpha))
        t_func_main[i] = mp.quad(lambda y: y * mp.exp(first_exp_term*y**alpha) * mp.exp(second_exp_term*y**2) * mp.exp(third_exp_term_main*y**2),[0,mp.inf])
        t_func_sc[i] = mp.quad(lambda y: y* mp.exp(first_exp_term*y**alpha) * mp.exp(second_exp_term*y**2) * mp.exp(third_exp_term_sc*y**2), [0,mp.inf])
    temp_second_sum = sum(2*mp.pi*density_update*t_func_sc)
    temp_third_sum = sum(2*mp.pi*density_org[0:-1]*t_func_main[0:-1])
    rate_coverage = (2*mp.pi*density_org[-1]/expected_activity)*(t_func_main[-1]*first_int + temp_second_sum*second_int) + temp_third_sum
    return (rate_coverage)

########################################################################################################################
def rate_cov_strate_sleeping_v3_uniform(k_matrix,alpha,rate_th1,lamda_u1,bw):
    # define the distribution of the activity
    first_int = (1/3)                                     # correspond to the integral in first term
    second_int = (1/2) - (1/3)                                            # correspond to the integral in second term
    #---------------------------- calibrated -----------------------------------------
    expected_activity = (1/2)                             # expected value of a      # this changes for optimization hance should be calibratable
    expected_strategic_function = (1/2)                           # expected value of s
    #---------------------------- calibrated -----------------------------------------
    noise_var = 1
    # preprocessing - define empty elements and get information about the inputs
    k_mat = k_matrix.copy()
    num_tiers = k_mat.shape[0]
    density_org = k_mat[:,2].copy()
    power = gb.db2pow(k_mat[:,1])
    density_update = np.array(density_org*([1]*(num_tiers-1)+[expected_strategic_function]))
    # define necessary values
    area_org = np.zeros(num_tiers,float)                      # original association probability
    area_sc_update = np.zeros(num_tiers,float)                # association probability of disconnected cell
    N_k_u1 = np.zeros(num_tiers,float)                        #number of users in tier K BS
    N_k_sc = np.zeros(num_tiers,float)                        #number of users in tier K BS
    N_k_total = np.zeros(num_tiers,float)
    threshold_u1 = np.zeros(num_tiers,float)                  #threshold for users in tier K BS
    t_func_main = np.zeros(num_tiers,float)
    t_func_sc = np.zeros(num_tiers,float)
    for i in range(num_tiers):
        area_org[i] = A_k(density_org,power,alpha,i)
        area_sc_update[i] = A_k(density_update,power,alpha,i)
    N_k_u1 = 1 + 1.28*lamda_u1*(area_org/density_org)
    N_k_sc = expected_activity*(1-expected_strategic_function)*density_org[-1]*N_k_u1[-1]*(area_sc_update/density_update)   #for binary optimization
    N_k_total = N_k_u1 + N_k_sc
    threshold_u1 = 2**((rate_th1/bw)*N_k_total) -1
    for i in range(num_tiers):
        first_exp_term = -(threshold_u1[i]*noise_var/power[i])
        z_term = 0 if (threshold_u1[i]==0) else (threshold_u1[i])**(2/alpha) *mp.quad(lambda u: 1/ (1+u**(alpha/2)),[(1/threshold_u1[i])**(2/alpha),mp.inf])
        second_exp_term = -mp.pi*z_term* sum(density_update*(power/power[i])**(2/alpha))
        third_exp_term_main = -mp.pi * sum(density_org*(power/power[i])**(2/alpha))
        third_exp_term_sc = -mp.pi * sum(density_update*(power/power[i])**(2/alpha))
        t_func_main[i] = mp.quad(lambda y: y * mp.exp(first_exp_term*y**alpha) * mp.exp(second_exp_term*y**2) * mp.exp(third_exp_term_main*y**2),[0,mp.inf])
        t_func_sc[i] = mp.quad(lambda y: y* mp.exp(first_exp_term*y**alpha) * mp.exp(second_exp_term*y**2) * mp.exp(third_exp_term_sc*y**2), [0,mp.inf])
    temp_second_sum = sum(2*mp.pi*density_update*t_func_sc)
    temp_third_sum = sum(2*mp.pi*density_org[0:-1]*t_func_main[0:-1])
    #rate_coverage = (2*mp.pi*density_org[-1]/expected_activity)*(t_func_main[-1]*first_int + temp_second_sum*second_int) + temp_third_sum
    rate_coverage = (area_org[-1]/expected_activity)*((2*mp.pi*density_org[-1]/area_org[-1])*t_func_main[-1]*first_int + temp_second_sum*second_int) + temp_third_sum
    #print((sum(2*mp.pi*density_update*t_func_sc)*second_int+t_func_main[-1]*first_int)*(2*mp.pi*density_org[-1]))
    #print(t_func_main[-1]*first_int)
    #print(temp_second_sum*second_int)
    #print((2*mp.pi*density_org[-1]/expected_activity)*(t_func_main[-1]*first_int) + temp_third_sum)
    #print(temp_second_sum)
    return (rate_coverage)
########################################################################################################################
'''
alpha = 4
lamda_u = 50
bw = 10*10**6
rate_th = 1 * 10**6
result = rate_cov_strate_sleeping_v3_uniform(k_mat1,alpha,rate_th,lamda_u,bw)
print(result)'''