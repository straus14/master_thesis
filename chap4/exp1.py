from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
import matplotlib.pyplot as plt
import random_sleeping as rs
import strategic_sleeping as ss
import pylab as pl
import opti_density_random_sleeping as opt1
########################################################################################################################

########################################################################################################################
            # define the system
########################################################################################################################
#-----------------------------------------------------------------------------------------------------------------------
#                 0            1             2               3
#                [K-tier,    Pow(dBm),     density,         BW ]
k1 =    np.array([1,         43,             1,          10*10**6],float)
k2 =    np.array([2,         38,             5,          10*10**6],float)
k3 =    np.array([3,         21,             10,         10*10**6],float)
k_mat1 = np.array([k1, k2, k3],float)

bw = 10*10**6
del1 = 4.7
del2 = 2.6
del3 = 4
b1_static = 130
b2_static = 56
b3_static = 6.8
b1_tx = gb.dbm2pow(43)
b2_tx = gb.dbm2pow(38)
b3_tx = gb.dbm2pow(21)
b1_ec = b1_tx + del1*b1_static*bw       # energy consumption at one BS in tier1
b2_ec = b2_tx + del2*b2_static*bw       # energy consumption at one BS in tier2
b3_ec = b3_tx + del3*b3_static*bw       # energy consumption at one BS in tier3
########################################################################################################################

'''
########################################################################################################################
          #chap32 - figure1 - rate coverage increases monotonically with BS density
########################################################################################################################
alpha = 4
rate_th1 = 1*10**6
rate_th2 = 0.5*10**6
lamda_u1 = 30
lamda_u2 = 20
bw = 10*10**6
q_ary = pl.frange(0,1,0.1)
rate_cov1 = np.zeros_like(q_ary,float)
rate_cov2 = np.zeros_like(q_ary,float)
rate_cov = np.zeros_like(q_ary,float)
kappa = (rate_th1*lamda_u1)/(rate_th1*lamda_u1 + rate_th2*lamda_u2)
for idx,val in enumerate(q_ary):
    rate_cov1[idx] = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th1,lamda_u1,val,kappa*bw)
    rate_cov1[idx] = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th2,lamda_u2,val,(1-kappa)*bw)
rate_cov = (lamda_u1/(lamda_u1+lamda_u2))*rate_cov1 + (lamda_u2/(lamda_u1+lamda_u2))*rate_cov2
plt.plot(q_ary,rate_cov,'b*-')
plt.xlabel('SC activation probability ($q$)')
plt.ylabel('rate coverage - $\mathcal{R}_{rs}(q)$')
plt.show()
#print(rate_cov)
########################################################################################################################
'''
'''
########################################################################################################################
         # --------------------- DuMMy-------------------NegleCt-----------------------------------------
########################################################################################################################
          #chap32 - figure2 - Energy Efficiency vs Small Cell Density
########################################################################################################################
alpha = 4
rate_th1 = 0.5*10**6
rate_th2_ary = pl.frange(1.5,3)*10**6
rate_th2 = 2*10**6
lamda_u1 = 30
lamda_u2 = 20
bw = 10*10**6
kappa = (rate_th1*lamda_u1)/(rate_th1*lamda_u1 + rate_th2*lamda_u2)
p1 = (lamda_u1)/(lamda_u1 + lamda_u2)
lamda_sc = np.arange(0,15,1)
rate_cov1_rs = np.zeros_like(lamda_sc,float)
rate_cov2_rs = np.zeros_like(lamda_sc,float)
EE_consum_rs = np.zeros_like(lamda_sc,float)
for idx,val in enumerate(lamda_sc):
    k_mat1[2][2] = val
    rate_cov1_rs[idx] = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th1,lamda_u1,1,kappa*bw)
    rate_cov2_rs[idx] = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th2,lamda_u2,1,(1-kappa)*bw)
    EE_consum_rs[idx] = b1_ec* k_mat1[0][2] + b2_ec* k_mat1[1][2] + b3_ec* val
ASE = ((rate_th1*lamda_u1*rate_cov1_rs)+(rate_th2*lamda_u2*rate_cov2_rs))
EE = ASE /EE_consum_rs
plt.plot(lamda_sc,EE,'b*-')
plt.xlabel('SCBS density')
plt.ylabel('Energy Efficiecy')
plt.show()
'''
'''
########################################################################################################################
          #chap32 - figure2 - Energy Efficiency vs Small Cell Density for various rate threshold of community 2
########################################################################################################################
alpha = 4
rate_th1 = 0.5*10**6
rate_th2_ary = pl.frange(1.5,3,0.5)*10**6

lamda_u1 = 30
lamda_u2 = 20
bw = 10*10**6
p1 = (lamda_u1)/(lamda_u1 + lamda_u2)
lamda_sc = np.arange(0,40,2)
EE_mat = np.zeros((len(rate_th2_ary),len(lamda_sc)))

for idx,val in enumerate(rate_th2_ary):
    rate_th2 = val
    kappa = (rate_th1*lamda_u1)/(rate_th1*lamda_u1 + rate_th2*lamda_u2)
    rate_cov1_rs = np.zeros_like(lamda_sc,float)
    rate_cov2_rs = np.zeros_like(lamda_sc,float)
    EE_consum_rs = np.zeros_like(lamda_sc,float)
    for idx2,val2 in enumerate(lamda_sc):
        k_mat1[2][2] = val2
        rate_cov1_rs[idx2] = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th1,lamda_u1,1,kappa*bw)
        rate_cov2_rs[idx2] = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th2,lamda_u2,1,(1-kappa)*bw)
        EE_consum_rs[idx2] = b1_ec* k_mat1[0][2] + b2_ec* k_mat1[1][2] + b3_ec* val2
    ASE = ((rate_th1*lamda_u1*rate_cov1_rs)+(rate_th2*lamda_u2*rate_cov2_rs))
    EE_mat[idx] = ASE /EE_consum_rs
np.savez('/home/straus14/hiwi_python/task_1/sleep_strategies/fig2',lamda_sc=lamda_sc,rate_th2_ary=rate_th2_ary,EE_mat=EE_mat)
#print(EE_mat)
########################################################################################################################
'''
'''
########################################################################################################################
          #chap32 - figure3 - variation of \lambda_{sc} vs \rho_{2}
########################################################################################################################
k12 = np.array([1,         43,             1,          10*10**6],float)
k22 = np.array([2,         38,             5,          10*10**6],float)
k32 = np.array([3,         21,             1,          10*10**6],float)
k_mat2 = np.array([k12, k22, k32],float)
alpha = 4

rate_th1 = 1 * 10**6
rate_th2_aray = pl.frange(1,5,1)*10**6
lamda_u1 = 25
lamda_u2 = 25
tar_cov1 = 0.3
bw = 10*10**6
tar_cov2_ary = np.array([0.5,0.4,0.3])
result_lamda_sc = np.zeros((len(tar_cov2_ary),len(rate_th2_aray)))
for idx1,val1 in enumerate(tar_cov2_ary):
    tar_cov2 = val1
    for idx2, val2 in enumerate(rate_th2_aray):
        rate_th2 = val2
        result_lamda_sc[idx1][idx2] = opt1.q_opt_rate_cov_fairness_v14(k_mat2,alpha,rate_th1,rate_th2,lamda_u1,lamda_u2,bw,tar_cov1,tar_cov2)
np.savez('/home/straus14/hiwi_python/task_1/sleep_strategies/fig3',tar_cov2_ary=tar_cov2_ary,rate_th2_aray=rate_th2_aray,result_lamda_sc=result_lamda_sc)
########################################################################################################################
'''
'''
########################################################################################################################
          #chap32 - figure4 - variation of \lambda_{sc} vs \rho_{2}
########################################################################################################################
k12 = np.array([1,         43,             1,          10*10**6],float)
k22 = np.array([2,         38,             5,          10*10**6],float)
k32 = np.array([3,         21,             1,          10*10**6],float)
k_mat2 = np.array([k12, k22, k32],float)
alpha = 4

rate_th1 = 1 * 10**6
rate_th2_aray = pl.frange(2,4,1)*10**6
lamda_u1 = 25
lamda_u2 = 25
tar_cov1 = 0.3
bw = 10*10**6
#tar_cov2_ary = np.array([0.2,0.7])
tar_cov2_ary = pl.frange(0.2,0.7,0.1)
result_lamda_sc = np.zeros((len(rate_th2_aray),len(tar_cov2_ary)))
for idx1, val1 in enumerate(rate_th2_aray):
    rate_th2 = val1
    for idx2,val2 in enumerate(tar_cov2_ary):
        tar_cov2 = val2
        result_lamda_sc[idx1][idx2] = opt1.q_opt_rate_cov_fairness_v14(k_mat2,alpha,rate_th1,rate_th2,lamda_u1,lamda_u2,bw,tar_cov1,tar_cov2)
np.savez('/home/straus14/hiwi_python/task_1/sleep_strategies/fig4',tar_cov2_ary=tar_cov2_ary,rate_th2_aray=rate_th2_aray,result_lamda_sc=result_lamda_sc)
########################################################################################################################
'''
'''
########################################################################################################################
          #chap32 - figure5 - variation of \lambda_{sc} vs \rho_{2}
########################################################################################################################
k12 = np.array([1,         43,             1,          10*10**6],float)
k22 = np.array([2,         38,             5,          10*10**6],float)
k32 = np.array([3,         21,             1,          10*10**6],float)
k_mat2 = np.array([k12, k22, k32],float)
alpha = 4

rate_th1 = 1 * 10**6
rate_th2_aray = pl.frange(2,4,1)*10**6
lamda_u1 = 25
lamda_u2_ary = pl.frange(25,55,10)
tar_cov1 = 0.5
tar_cov2 = 0.5
bw = 10*10**6
result_lamda_sc = np.zeros((len(rate_th2_aray),len(lamda_u2_ary)))
for idx1, val1 in enumerate(rate_th2_aray):
    rate_th2 = val1
    for idx2, val2 in enumerate(lamda_u2_ary):
        lamda_u2 = val2
        result_lamda_sc[idx1][idx2] = opt1.q_opt_rate_cov_fairness_v14(k_mat2,alpha,rate_th1,rate_th2,lamda_u1,lamda_u2,bw,tar_cov1,tar_cov2)
np.savez('/home/straus14/hiwi_python/task_1/sleep_strategies/fig5',lamda_u2_ary=lamda_u2_ary,rate_th2_aray=rate_th2_aray,result_lamda_sc=result_lamda_sc)
########################################################################################################################
'''
'''
########################################################################################################################
          #chap32 - figure6 - variation of \lambda_{sc} vs \rho_{2}
########################################################################################################################
k12 = np.array([1,         43,             1,          10*10**6],float)
k22 = np.array([2,         38,             5,          10*10**6],float)
k32 = np.array([3,         21,             1,          10*10**6],float)
k_mat2 = np.array([k12, k22, k32],float)
alpha = 4

rate_th1 = 1 * 10**6
rate_th2 = 2* 10**6
tar_cov1 = 0.3
tar_cov2_ary = np.array([0.3,0.5,0.7])
bw = 10*10**6
lamda_u1 = 25
lamda_u2_ary = pl.frange(25,50,5)
result_lamda_sc = np.zeros((len(tar_cov2_ary),len(lamda_u2_ary)))
for idx1,val1 in enumerate(tar_cov2_ary):
    tar_cov2 = val1
    for idx2, val2 in enumerate(lamda_u2_ary):
        lamda_u2 = val2
        result_lamda_sc[idx1][idx2] = opt1.q_opt_rate_cov_fairness_v14(k_mat2,alpha,rate_th1,rate_th2,lamda_u1,lamda_u2,bw,tar_cov1,tar_cov2)
np.savez('/home/straus14/hiwi_python/task_1/sleep_strategies/fig6',lamda_u2_ary=lamda_u2_ary,tar_cov2_ary=tar_cov2_ary,result_lamda_sc=result_lamda_sc)
########################################################################################################################
'''
'''
########################################################################################################################
          #chap32 - figure7 - venergy efficiency comparision of SS and RS
########################################################################################################################
alpha = 4
rate_th1 = 2 * 10**6
rate_th2 = 4 * 10**6
lamda_u1 = 25
lamda_u2 = 25
bw = 10*10**6
kappa = (lamda_u1*rate_th1)/(lamda_u1*rate_th1 + lamda_u2*rate_th2)
p1 = (lamda_u1)/(lamda_u1+lamda_u2)
q_on = pl.frange(0.1,1,0.1)
energy_consum = np.zeros(len(q_on))
rate_cov_rs1 = np.zeros(len(q_on))
rate_cov_rs2 = np.zeros(len(q_on))
rate_cov_ss1 = np.zeros(len(q_on))
rate_cov_ss2 = np.zeros(len(q_on))
ee_total = np.zeros(len(q_on))
for idx1, val1 in enumerate(q_on):
    rate_cov_rs1[idx1] = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th1,lamda_u1,val1,kappa*bw)
    rate_cov_rs2[idx1] = rs.rate_cov_random_switching_v2(k_mat1,alpha,rate_th2,lamda_u2,val1,(1-kappa)*bw)
    rate_cov_ss1[idx1] = ss.rate_cov_strate_sleeping_v2(k_mat1,alpha,rate_th1,lamda_u1,kappa*bw,val1)
    rate_cov_ss2[idx1] = ss.rate_cov_strate_sleeping_v2(k_mat1,alpha,rate_th2,lamda_u2,(1-kappa)*bw,val1)
    ee_total[idx1] = b1_ec*k_mat1[0][2] + b2_ec*k_mat1[1][2] + b3_ec*k_mat1[0][2]*val1 + (1-val1)*4.3
ee_rs = (lamda_u1*rate_th1*rate_cov_rs1 + lamda_u2*rate_th2*rate_cov_rs2)/ee_total
ee_ss = (lamda_u1*rate_th1*rate_cov_ss1 + lamda_u2*rate_th2*rate_cov_ss2)/ee_total
np.savez('/home/straus14/hiwi_python/task_1/sleep_strategies/fig7',q_on=q_on,ee_rs=ee_rs,ee_ss=ee_ss)
########################################################################################################################
'''
'''
########################################################################################################################
          #chap32 - figure8 - rate coverage SS - binary vs uniform activity distribution - single community user
########################################################################################################################
alpha = 4
lamda_u = 50
bw = 10*10**6
rate_th_ary = pl.frange(0,4,0.5)*10**6
rate_cov_ss_binary = np.zeros(len(rate_th_ary))
rate_cov_ss_uniform = np.zeros(len(rate_th_ary))
bin_q = 0.5
for idx,val in enumerate(rate_th_ary):
    rate_cov_ss_binary[idx] = ss.rate_cov_strate_sleeping_v3_bin(k_mat1,alpha,val,lamda_u,bw,bin_q)
    rate_cov_ss_uniform[idx] = ss.rate_cov_strate_sleeping_v3_uniform(k_mat1,alpha,val,lamda_u,bw)
#print(rate_th_ary)
#print(rate_cov_ss_binary)
np.savez('/home/straus14/hiwi_python/task_1/sleep_strategies/fig8',rate_th_ary=rate_th_ary,rate_cov_ss_binary=rate_cov_ss_binary,rate_cov_ss_uniform=rate_cov_ss_uniform)
########################################################################################################################
'''






########################################################################################################################