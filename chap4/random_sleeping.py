#random switching
from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
from scipy.integrate import dblquad,quad
import matplotlib.pyplot as plt
#######################################################################################################################

#area calculation for tiers with equal pathloss exponent
def A_k(density,power,alpha,idx):
    area_k = density[idx] / (sum(density*(power/power[idx])**(2/alpha)))
    return area_k

#-----------------------------------------------------------------------------------------------------------------------
#this function gives the rate coverage probability of K-tier system

def rate_cov_random_switching(k_matrix,alpha,rate_th,lamda_user):
    k_mat = k_matrix.copy()
    num_tiers = k_mat.shape[0]
    density = k_mat[:,2]
    power = gb.db2pow(k_mat[:,1])
    bandwidth = k_mat[:,3]
    small_cell_idx = num_tiers -1                                                 #indicates the index of the small cell
    # initialize the integration result matrix
    tier_integ_result = np.zeros(num_tiers,float)                                      #integration results of each tier
    area_tiers = np.zeros(num_tiers,float)                                                            #area of the tiers
    threshold_tier = np.zeros(num_tiers,float)                                                   #threshold of the tiers
    N_k = np.zeros(num_tiers,float)                                                        #number of users in each tier
    first_exp_term = np.zeros(num_tiers,float)                                       #first exponential term in integral
    second_exp_term = np.zeros(num_tiers,float)                                     #second exponential term in integral
    third_exp_term = np.zeros(num_tiers,float)                                       #third exponential term in integral
    for i in range(num_tiers):
        area_tiers[i] = A_k(density,power,alpha,i)
        N_k[i] = 0 if density[i]==0 else 1.28*lamda_user*area_tiers[i]/density[i]
        threshold_tier[i] = 2**(rate_th*N_k[i]/bandwidth[i]) - 1
        first_exp_term[i] = threshold_tier[i]*1/power[i]
        third_exp_term[i] = mp.pi*sum(density*(power/power[i])**(2/alpha))
        Z_term = 0 if (threshold_tier[i]==0) else threshold_tier[i]**(2/alpha) * mp.quad(lambda u: 1 / (1 + u**(alpha/2)), [(1/threshold_tier[i])**(2/alpha), mp.inf])
        second_exp_term[i] = third_exp_term[i] * Z_term
    for k in range(num_tiers):
        tier_integ_result[k] = mp.quad(lambda y: y * mp.exp(-first_exp_term[k]*y**alpha) * mp.exp(-second_exp_term[k]*y**2) * mp.exp(-third_exp_term[k]* y**2), [0, mp.inf])
    rate_cov_prob = 2*mp.pi*sum(density*tier_integ_result)
    #print(rate_cov_prob)
    return (rate_cov_prob)
#-----------------------------------------------------------------------------------------------------------------------
#                 0            1             2               3
#                [K-tier,    Pow(dBm),     density,         BW ]
k1 =    np.array([1,         53,             1,          10*10**6],float)
k2 =    np.array([2,         33,             5,          10*10**6],float)
k3 =    np.array([3,         23,             10,         10*10**6],float)
k_mat1 = np.array([k1, k2, k3],float)
alpha = 4
rate_th = 1 * 10**6
lamda_user = 50
#rate_cov = rate_cov_random_switching(k_mat1,alpha,rate_th,lamda_user)
#print(rate_cov)
#executing and plotting result of rate coverage probability for 3-tier system
#-----------------------------------------------------------------------------------------------------------------------
'''rate_th1 = (np.array(range(0,21,2),float)/10)*10**6
rate_cov = np.zeros_like(rate_th1,float)
for idx, val in enumerate(rate_th1):
    rate_cov[idx] = rate_cov_random_switching(k_mat1,alpha,val,lamda_user)
#print(rate_cov)
#plotting of the result
plt.plot(rate_th1/10**6,rate_cov,'b*-')
plt.xticks(np.arange(0,2.2,0.2))
plt.title("rate coverage probabilty K-tier system")
plt.xlabel('rate threshold(Mbps)')
plt.ylabel('rate coverage probability')
plt.show()'''
#-----------------------------------------------------------------------------------------------------------------------

#plot for rate coverage always increases with density of small cell base station
#-----------------------------------------------------------------------------------------------------------------------
'''q_on = np.array(range(0,11),float)/10
rate_coverage = np.zeros_like(q_on,float)
for idx, val in enumerate(q_on):
    k_mat1[2][2] = val * k3[2]
    rate_coverage[idx] = rate_cov_random_switching(k_mat1,alpha,rate_th,lamda_user)
#plotting of the result
plt.plot(q_on,rate_coverage,'b*-')
plt.xticks(np.arange(0,1,0.1))
plt.title("rate coverage increases monotonically with density of small cell activation")
plt.xlabel('activation probability(q_on)')
plt.ylabel('rate coverage probability')
plt.show()'''



#-----------------------------------------------------------------------------------------------------------------------
#updated function with probability and bandwidth of small cell
# this function gives rate coverage when we specify small cell active probability 'q' and bandwidth 'bw'
# note small cell is considered to give at the end of the k-tiers in input
# this is for single community
#-----------------------------------------------------------------------------------------------------------------------

def rate_cov_random_switching_v2(k_matrix,alpha,rate_th,lamda_user,q_on,bw):
    k_mat = k_matrix.copy()
    num_tiers = k_mat.shape[0]
    #density = k_mat[:,2]
    #density[-1] = q_on
    density = np.array(k_mat[:,2]*([1]*(num_tiers-1)+[q_on]))                    #hence last row always corresponds to small cells
    #print(density)
    power = gb.db2pow(k_mat[:,1])
    small_cell_idx = num_tiers -1                                                 #indicates the index of the small cell
    #density[small_cell_idx] = q_on
    # initialize the integration result matrix
    tier_integ_result = np.zeros(num_tiers,float)                                      #integration results of each tier
    area_tiers = np.zeros(num_tiers,float)                                                            #area of the tiers
    threshold_tier = np.zeros(num_tiers,float)                                                   #threshold of the tiers
    N_k = np.zeros(num_tiers,float)                                                        #number of users in each tier
    first_exp_term = np.zeros(num_tiers,float)                                       #first exponential term in integral
    second_exp_term = np.zeros(num_tiers,float)                                     #second exponential term in integral
    third_exp_term = np.zeros(num_tiers,float)                                       #third exponential term in integral
    for i in range(num_tiers):
        area_tiers[i] = A_k(density,power,alpha,i)
        #N_k[i] = 0 if density[i]==0 else  1.28*lamda_user*area_tiers[i]/density[i]
        N_k[i] = 0 if density[i]==0 else 1 + 1.28*lamda_user*area_tiers[i]/density[i]
        threshold_tier[i] = mp.inf if (bw==0) else 2**(rate_th*N_k[i]/bw) - 1
        first_exp_term[i] = threshold_tier[i]*1/power[i]
        third_exp_term[i] = mp.pi*sum(density*(power/power[i])**(2/alpha))
        Z_term = 0 if (threshold_tier[i]==0) else threshold_tier[i]**(2/alpha) * mp.quad(lambda u: 1 / (1 + u**(alpha/2)), [(threshold_tier[i])**(-2/alpha), mp.inf])
        second_exp_term[i] = third_exp_term[i] * Z_term
    for k in range(num_tiers):
        tier_integ_result[k] = mp.quad(lambda y: y * mp.exp(-first_exp_term[k]*y**alpha) * mp.exp(-second_exp_term[k]*y**2) * mp.exp(-third_exp_term[k]* y**2), [0, mp.inf])
    rate_cov_prob = 2*mp.pi*sum(density*tier_integ_result)
    return (rate_cov_prob)

#-----------------------------------------------------------------------------------------------------------------------
                                                   #testing
#-----------------------------------------------------------------------------------------------------------------------
'''bw = 10 * 10**6
rate_threshold = np.array(range(0,21))/10 * 10**6
rate_cov = np.zeros_like(rate_threshold,float)
for idx, val in enumerate(rate_threshold):
    rate_cov[idx] = rate_cov_random_switching_v2(k_mat1,alpha,val,lamda_user,1,bw)
plt.plot(rate_threshold,rate_cov,'b*-')
plt.show()'''


'''
l_u1 = 25
l_u2 = 25
r1 = 4*10**6
r2 = 2*10**6
p1 = l_u1 / (l_u1 + l_u2)
p_c1 = (l_u1*r1) / (l_u1*r1 + l_u2*r2)
bw = 10 * 10**6
rate_cov1 = rate_cov_random_switching_v2(k_mat1,alpha,r1,l_u1,1,p_c1*bw)
rate_cov2 = rate_cov_random_switching_v2(k_mat1,alpha,r2,l_u2,1,(1-p_c1)*bw)
rate_cov_multiple = p1 * rate_cov1 + (1-p1)*rate_cov2
rate_cov_single = rate_cov_random_switching_v2(k_mat1,alpha,max(r1,r2),(l_u1+l_u2),1,bw)
print(rate_cov_single)
print(rate_cov_multiple)
'''

#-----------------------------------------------------------------------------------------------------------------------
#                 0            1             2               3
#                [K-tier,    Pow(dBm),     density,         BW ]
k12 =    np.array([1,         43,             1,          10*10**6],float)
k22 =    np.array([2,         38,             5,          10*10**6],float)
k32 =    np.array([3,         21,             10,         10*10**6],float)
k_mat2 = np.array([k12, k22, k32],float)
alpha = 4
rate_th1 = 0.5 * 10**6
rate_th2 = 1 * 10**6
lamda_u1 = 25
lamda_u2 = 25
#lamda_user = 50
bw = 10*10**6
q_on = 1
'''
result = rate_cov_random_switching_v2(k_mat2,alpha,rate_th1,lamda_u1,q_on,bw)
print(result)'''