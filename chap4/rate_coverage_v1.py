
#specific functions for the IEEE paper:"Offloading in Heterogeneous Network: Modeling, Analysis, and Design Insights"
from __future__ import division
import numpy as np
import mpmath as mp
import general_library as gb
import matplotlib.pyplot as plt
##############################################################
       #to check for execution time (debugging)
import time
##############################################################
 
#G_ij(m,k) function ------------equation-7
def G_mk(density,assoc_weight_norm,alpha):
    density_array = np.array(density)
    assoc_weight_norm_array = np.array(assoc_weight_norm)
    alpha_array = np.array(alpha)
    return (density_array*(assoc_weight_norm_array**(2/alpha_array)))
 
#A_ij--- probability a user associated with an AP ---------- equation-6
def A_ij(density,G_ij,alpha_norm):
    f1 = lambda z: sum((G_ij)*(z**(2/alpha_norm)))
    return (mp.quad(lambda z:(z * mp.exp(-mp.pi * f1(z))),[0,mp.inf])*(2*mp.pi*density))
 
# Z(a,b,c) function--------------Lemma-5
def Z_abc(a,b,c):
    if(a!=0):
        return(mp.quad(lambda u:(1/(1+u**(b/2))),[((c/a)**(2/b)),mp.inf])*(a**(2/b)))
    else:
        return (0)
 
#D_ij(k,T_ij) function----------Lemma-5 for theorem 1 and Corollary 1
#input format matrix for tiers elements of corresponding RAT system, where each row is       [P^_ik, T^_ik, alpha_k, L_k, L"_k]---(Kx5 matrix)
#and fixed parameter t(SINR threshold of the cell)
def D_ij(matrix_ik,t):
    out = np.zeros(matrix_ik.shape[0],float)#initialize the output array
    for i in range(len(out)):
        if(matrix_ik[i,0]):
            out[i] = (matrix_ik[i,0]**(2/matrix_ik[i,2]))*((matrix_ik[i,3]*Z_abc(t,matrix_ik[i,2],matrix_ik[i,1]/matrix_ik[i,0]))+(matrix_ik[i,4]*Z_abc(t,matrix_ik[i,2],0)))
    return(out)
 
#t_x function------subfunction in Theorem1
def t_x(x):
    return ((2**x)-1)
#                                                              0      1       2        3     4     5
#D_ij (k,T_ij) function----------Lemma-5 for SINR threshold  [P^_ik, T^_ik, alpha_k, t_ij,  L_k, L"_k]---(Kx5 matrix)
def Dij_normal(k_matrix):
    out = np.zeros(k_matrix.shape[0],float)#initialize the output array
    for idx,val in enumerate(out):
        if(k_matrix[idx,0]):
            out[idx] = (k_matrix[idx,4]*Z_abc(k_matrix[idx,3],k_matrix[idx,2],k_matrix[idx,1]/k_matrix[idx,0])+k_matrix[idx,5]*Z_abc(k_matrix[idx,3],k_matrix[idx,2],0))*(k_matrix[idx,0]**(2/k_matrix[idx,2]))
    return (out)
 
#-----------------------------------------------------------------------------------------------------------------------
                                    #calibration parameters
max_num_tiers_cl = 3
 
 
#-----------------------------------------------------------------------------------------------------------------------
 
#-----------------------------------------------------------------------------------------------------------------------
                                           #Corollary1-----load at each AP is equal---
#-----------------------------------------------------------------------------------------------------------------------
#Inputs for the rate coverage probability
#[M-rat, K-tier, openaccess_flag, Pow(dBm), density, alpha, Bias(dB), BW ]----one row of (mk x 8) matrix
#(density of the users)
# rate threshold for which probability is needed
def rate_coverage_meanload(mk_matrix,lamda_u,rate_threshold,bw):
    mk_mat = mk_matrix.copy()
    mk_mat_size = mk_mat.shape
    bw = 10*10**6 if (bw==0) else bw
    mk_mat[:,3] = gb.dbm2pow(mk_mat[:,3])     # converts power from dBm to W
    mk_mat[:,6] = gb.db2pow(mk_mat[:,6])      # coverts Bias from dB to number
    mk_OpenCell_RowIdx = (np.asarray(np.where(mk_mat[:,2]==True))).flatten() # gives row index of open-access (m,k)
    area_opencells = np.zeros(mk_mat_size[0],float)  # initialize area of open cells to zero
    Nij = np.zeros(mk_mat_size[0],float)  #initialize the E[N_ij] value to zero
    Gij_mk_array = [[] for i in range(mk_mat_size[0])]                            #Gij(m,k) polynomial coefficient matrix
    alpha_norm_mk_array = [[] for i in range(mk_mat_size[0])]
    Tmk = mk_mat[:,3]*mk_mat[:,6]                #Association weight
    Tmk_OpenCell = Tmk[mk_OpenCell_RowIdx]        #Association weight for open access cells
    temp_integral_result = np.zeros(mk_mat_size[0],float)
    #rate_coverage = np.zeros_like(rate_threshold,float)   #output array
    for i in mk_OpenCell_RowIdx:            # loop to calculate Area, Nij and Gij(m,k) of open-access cells
        Tmk_norm = Tmk_OpenCell/Tmk[i]
        alpha_norm = mk_mat[mk_OpenCell_RowIdx,5]/mk_mat[i,5]
        Gij_mk = mk_mat[mk_OpenCell_RowIdx,4]*(Tmk_norm**(2/mk_mat[mk_OpenCell_RowIdx,5]))  # Eq-7 is computed with (i,j) corresponding to ith row index
        Gij_mk_array[i] = Gij_mk
        alpha_norm_mk_array[i] = alpha_norm
        area_opencells[i] = A_ij(mk_mat[i,4],Gij_mk,alpha_norm)
        Nij[i] = 0 if (mk_mat[i,4]==0) else (1.28*lamda_u*area_opencells[i]/mk_mat[i,4])+1
    #for rate_idx in range(len(rate_threshold)):
        for i in mk_OpenCell_RowIdx:
            Dij_inp = np.array([np.zeros(5,float) for tier in range(3)])    #initialize matrix for Dij(m,k) function
            tier_RowIndx = (np.asarray(np.where(mk_mat[:,0]==mk_mat[i,0]))).flatten() #index of tiers of ith rows RAT system
            for x in tier_RowIndx:                    #fill Dij(m,k) input matrix
                k = mk_mat[x,1]                       # gives the Kth tier of Mth RAT system
                Dij_inp[k-1,0]=mk_mat[x,3]/mk_mat[i,3]
                Dij_inp[k-1,1]=Tmk[x]/Tmk[i]
                Dij_inp[k-1,2]=mk_mat[x,5]
                if(mk_mat[x,2]):
                    Dij_inp[k-1,3]=mk_mat[x,4]
                else:
                    Dij_inp[k-1,4]=mk_mat[x,4]
            Dij_poly = D_ij(Dij_inp,t_x((rate_threshold/bw)*Nij[i]))
            alpha_norm_dij = Dij_inp[:,2]/mk_mat[i,5]
            alpha_norm_dij[alpha_norm_dij==0]=1  #to avoid divide by zero warning :):):)
            f1 = lambda y: sum(Dij_poly*(y**(2/alpha_norm_dij)))
            f2 = lambda y: sum(np.asarray(Gij_mk_array[i])*(y**(2/np.asarray(alpha_norm_mk_array[i]))))
            f3 = lambda y: (t_x((rate_threshold/bw)*Nij[i])*y**(mk_mat[i,5]))/mk_mat[i,3]
            f4 = lambda y: -f3(y)-mp.pi*f2(y)-mp.pi*f1(y)
            temp_integral_result[i] = mp.quad(lambda y: y * mp.exp(f4(y)),[0,mp.inf])
        rate_coverage=(sum((2*np.pi*mk_mat[:,4])*temp_integral_result))
    return (rate_coverage)
 
#-----------------------------------------------------------------------------------------------------------------------
                                                #Theorem_1
#-----------------------------------------------------------------------------------------------------------------------
 
 
def rate_coverage(mk_matrix,lamda_u,rate_threshold):
    mk_mat = mk_matrix.copy()
    mk_mat_size = mk_mat.shape
    mk_mat[:,3] = gb.dbm2pow(mk_mat[:,3])     # converts power from dBm to W
    mk_mat[:,6] = gb.db2pow(mk_mat[:,6])      # coverts Bias from dB to number
    mk_OpenCell_RowIdx = (np.asarray(np.where(mk_mat[:,2]==True))).flatten() # gives row index of open-access (m,k)
    area_opencells = np.zeros(mk_mat_size[0],float)  # initialize area of open cells to zero
    Gij_mk_array = [[] for i in range(mk_mat_size[0])]  #Gij(m,k) polynomial coefficient matrix
    alpha_norm_mk_array = [[] for i in range(mk_mat_size[0])]
    Tmk = mk_mat[:,3]*mk_mat[:,6]                 #Association weight
    Tmk_OpenCell = Tmk[mk_OpenCell_RowIdx]        #Association weight for open access cells
    temp_sum_integral_result = np.zeros(mk_mat_size[0],float)
    rate_coverage = np.zeros_like(rate_threshold,float)   #output array
    N_max = 4*lamda_u
    for i in mk_OpenCell_RowIdx:            # loop to calculate Area, Nij and Gij(m,k) of open-access cells
        Tmk_norm = Tmk_OpenCell/Tmk[i]
        alpha_norm = mk_mat[mk_OpenCell_RowIdx,5]/mk_mat[i,5]
        Gij_mk = mk_mat[mk_OpenCell_RowIdx,4]*(Tmk_norm**(2/mk_mat[mk_OpenCell_RowIdx,5]))  # Eq-7 is computed with (i,j) corresponding to ith row index
        Gij_mk_array[i] = Gij_mk
        alpha_norm_mk_array[i] = alpha_norm
        area_opencells[i] = A_ij(mk_mat[i,4],Gij_mk,alpha_norm)
    for rate_idx in range(len(rate_threshold)):
        for i in mk_OpenCell_RowIdx:
            f1_temp = lamda_u*area_opencells[i]/mk_mat[i,4]
            f1= lambda n:((3.5**3.5)/mp.fac(n))*(mp.gamma(n+4.5)/3.32335097044784)*(f1_temp**n)*(3.5+f1_temp)**(-n-4.5)  #computation of sum terms in equation 20
            temp_sum = np.zeros(N_max+1,float)
            temp_integral = np.zeros(N_max+1,float)
            f2 = lambda y: sum(np.asarray(Gij_mk_array[i])*(y**(2/np.asarray(alpha_norm_mk_array[i])))) #Gij(m,k) term in eq 20
            Dij_inp = np.array([np.zeros(5,float) for tier in range(3)])    #initialize matrix for Dij(m,k) function
            #alpha_norm_dij = Dij_inp[:,2]/mk_mat[i,5]
            #alpha_norm_dij[alpha_norm_dij==0]=1  #to avoid divide by zero warning :):):)
            tier_RowIndx = (np.asarray(np.where(mk_mat[:,0]==mk_mat[i,0]))).flatten() #index of tiers of ith rows RAT system
            for x in tier_RowIndx:                    #fill Dij(m,k) input matrix
                k = mk_mat[x,1]                       # gives the Kth tier of Mth RAT system
                Dij_inp[k-1,0]=mk_mat[x,3]/mk_mat[i,3]
                Dij_inp[k-1,1]=Tmk[x]/Tmk[i]
                Dij_inp[k-1,2]=mk_mat[x,5]
                if(mk_mat[x,2]):
                    Dij_inp[k-1,3]=mk_mat[x,4]
                else:
                    Dij_inp[k-1,4]=mk_mat[x,4]
            for n in range(0,N_max+1):
                temp_sum[n] = f1(n)
                sinr = (rate_threshold[rate_idx]/mk_mat[i,7])*(n+1)
                #print(sinr)
                Dij_poly = D_ij(Dij_inp,t_x((rate_threshold[rate_idx]/mk_mat[i,7])*(n+1)))
                alpha_norm_dij = Dij_inp[:,2]/mk_mat[i,5]
                alpha_norm_dij[alpha_norm_dij==0]=1  #to avoid divide by zero warning :):):)
                f3 = lambda y: sum(Dij_poly*(y**(2/alpha_norm_dij)))  #Dij term in integration
                f4 = lambda y: (t_x((rate_threshold[rate_idx]/mk_mat[i,7])*(n+1))*y**(mk_mat[i,5]))/mk_mat[i,3]
                f5 = lambda y: -f4(y)-mp.pi*f3(y)-mp.pi*f2(y)
                temp_integral[n]=mp.quad(lambda y: y * mp.exp(f5(y)),[0,mp.inf])
            temp_sum_integral_result[i] = np.dot(temp_sum,temp_integral)
        rate_coverage[rate_idx] = (sum((2*np.pi*mk_mat[:,4])*temp_sum_integral_result))
    return (rate_coverage,rate_threshold)
 
 
#-----------------------------------------------------------------------------------------------------------------------
                                  #examples below (just uncomment to execute particular set)
#-----------------------------------------------------------------------------------------------------------------------
 
 
#                     0       1            2              3        4          5          6          7
#                   [M-rat, K-tier, openaccess_flag, Pow(dBm), density,     alpha,    Bias(dB),     BW ]
m_1_k_1 =  np.array([1,       1,           True,          53,        1,         3.5,       0,        10*10**6],float)
m_2_k_3 =  np.array([2,       3,           True,          23,       10,         4,         15,        10*10**6],float)
m_2_k_3c = np.array([2,       3,           False,         23,       10,         4,         15,        10*10**6],float)
mk_mat1 = np.array([ m_1_k_1, m_2_k_3, m_2_k_3c])
mk_mat2 = mk_mat1.copy()
mk_mat2[:,6] = mk_mat2[:,6]*3
lamda_u = 50
rate_threshold1 = ((np.array(range(0,22,1),float)/10)*10**6)
rate_threshold2 = ((np.array(range(0,60,5),float)/10)*10**6)
rate_threshold3 = np.array([1*10**6])
 
#                     0       1            2              3        4          5          6          7
#                   [M-rat, K-tier, openaccess_flag, Pow(dBm), density,     alpha,    Bias(dB),     BW ]
m1k1 =  np.array([1,       1,           True,          53,        1,         3.5,       0,        10*10**6])
m1k2 =  np.array([1,       2,           True,          33,        5,         3.8,       5,        10*10**6])
m2k2 =  np.array([2,       2,           True,          33,        5,         3.8,       5,        10*10**6])
m2k3 =  np.array([2,       3,           True,          23,        10,         4,        10,       10*10**6])
mk_mat3 = np.array([m1k1,m1k2,m2k2,m2k3])
#mk_mat4 = mk_mat3.copy()
#mk_mat4[3,6] = 10
bw = 10*10**6

'''cr_p,cr_r = rate_coverage_meanload(mk_mat1,lamda_u,rate_threshold1,bw)
plt.plot(cr_r/10**6,cr_p,'r-*')
plt.show()'''

'''
th_p,th_r = rate_coverage(mk_mat1,lamda_u,rate_threshold1)
cr_p,cr_r = rate_coverage_meanload(mk_mat1,lamda_u,rate_threshold1)
plt.plot(th_r/10**6,th_p,'g-o',cr_r/10**6,cr_p,'r-^')
plt.xticks(np.arange(0,2.2,0.2))
plt.yticks(np.arange(0,1.1,0.1))
plt.title("rate probabilty vs rate threshold with figure2, B23 = 5dB")
plt.xlabel('rate threshold')
plt.ylabel('rate probability')
plt.legend(['theorem1','corollary1'])
plt.show()
print("rate probabilty vs rate threshold with figure2, B23 = 5dB")
print(th_r)
print(th_p)
print(cr_p)'''
 
'''th_p,th_r = rate_coverage(mk_mat2,lamda_u,rate_threshold1)
cr_p,cr_r = rate_coverage_meanload(mk_mat2,lamda_u,rate_threshold1)
plt.plot(th_r/10**6,th_p,'g-o',cr_r/10**6,cr_p,'r-^')
plt.xticks(np.arange(0,2.2,0.2))
plt.yticks(np.arange(0,1.1,0.1))
plt.title("rate probabilty vs rate threshold figure2, B23 = 15dB")
plt.xlabel('rate threshold')
plt.ylabel('rate probability')
plt.legend(['theorem1','corollary1'])
plt.show()
print("rate probabilty vs rate threshold figure2, B23 = 15dB")
print(th_r)
print(th_p)
print(cr_p)'''
 
'''th_p,th_r = rate_coverage(mk_mat3,lamda_u,rate_threshold2)
cr_p,cr_r = rate_coverage_meanload(mk_mat3,lamda_u,rate_threshold2)
plt.plot(th_r/10**6,th_p,'g-o',cr_r/10**6,cr_p,'r-^')
plt.xticks(np.arange(0,2.2,0.2))
plt.yticks(np.arange(0,1.1,0.1))
plt.title("rate probabilty vs rate threshold figure3, B23 = 0dB")
plt.xlabel('rate threshold')
plt.ylabel('rate probability')
plt.legend(['theorem1','corollary1'])
plt.show()
print("rate probabilty vs rate threshold figure3, B23 = 0dB")
print(th_r)
print(th_p)
print(cr_p)'''
 
'''th_p,th_r = rate_coverage(mk_mat3,lamda_u,rate_threshold2)
cr_p,cr_r = rate_coverage_meanload(mk_mat3,lamda_u,rate_threshold2)
plt.plot(th_r/10**6,th_p,'g-o',cr_r/10**6,cr_p,'r-^')
plt.xticks(np.arange(0,6,0.5))
plt.yticks(np.arange(0,1.1,0.1))
plt.title("rate probabilty vs rate threshold figure3, B23 = 10dB")
plt.xlabel('rate threshold')
plt.ylabel('rate probability')
plt.legend(['theorem1','corollary1'])
plt.show()
print("rate probabilty vs rate threshold figure3, B23 = 10dB")
print(th_r)
print(th_p)
print(cr_p)'''
 
#-----------------------------------------------------------------------------------------------------------------------
                                           #Corollary1-----load at each AP is equal---interference limited
#-----------------------------------------------------------------------------------------------------------------------
#Inputs for the rate coverage probability
#[M-rat, K-tier, openaccess_flag, Pow(dBm), density, alpha, Bias(dB), BW ]----one row of (mk x 8) matrix
#(density of the users)
# rate threshold for which probability is needed
def rate_coverage_meanload_int_lmt(mk_matrix,lamda_u,rate_threshold,bw):
    mk_mat = mk_matrix.copy()
    mk_mat_size = mk_mat.shape
    bw = 10*10**6 if (bw==0) else bw
    mk_mat[:,3] = gb.dbm2pow(mk_mat[:,3])     # converts power from dBm to W
    mk_mat[:,6] = gb.db2pow(mk_mat[:,6])      # coverts Bias from dB to number
    mk_OpenCell_RowIdx = (np.asarray(np.where(mk_mat[:,2]==True))).flatten() # gives row index of open-access (m,k)
    area_opencells = np.zeros(mk_mat_size[0],float)  # initialize area of open cells to zero
    Nij = np.zeros(mk_mat_size[0],float)  #initialize the E[N_ij] value to zero
    Gij_mk_array = [[] for i in range(mk_mat_size[0])]                            #Gij(m,k) polynomial coefficient matrix
    alpha_norm_mk_array = [[] for i in range(mk_mat_size[0])]
    Tmk = mk_mat[:,3]*mk_mat[:,6]                #Association weight
    Tmk_OpenCell = Tmk[mk_OpenCell_RowIdx]        #Association weight for open access cells
    temp_integral_result = np.zeros(mk_mat_size[0],float)
    #rate_coverage = np.zeros_like(rate_threshold,float)   #output array
    for i in mk_OpenCell_RowIdx:            # loop to calculate Area, Nij and Gij(m,k) of open-access cells
        Tmk_norm = Tmk_OpenCell/Tmk[i]
        alpha_norm = mk_mat[mk_OpenCell_RowIdx,5]/mk_mat[i,5]
        Gij_mk = mk_mat[mk_OpenCell_RowIdx,4]*(Tmk_norm**(2/mk_mat[mk_OpenCell_RowIdx,5]))  # Eq-7 is computed with (i,j) corresponding to ith row index
        Gij_mk_array[i] = Gij_mk
        alpha_norm_mk_array[i] = alpha_norm
        area_opencells[i] = A_ij(mk_mat[i,4],Gij_mk,alpha_norm)
        Nij[i] = 0 if (mk_mat[i,4]==0) else (1.28*lamda_u*area_opencells[i]/mk_mat[i,4])+1
    #for rate_idx in range(len(rate_threshold)):
        for i in mk_OpenCell_RowIdx:
            Dij_inp = np.array([np.zeros(5,float) for tier in range(3)])    #initialize matrix for Dij(m,k) function
            tier_RowIndx = (np.asarray(np.where(mk_mat[:,0]==mk_mat[i,0]))).flatten() #index of tiers of ith rows RAT system
            for x in tier_RowIndx:                    #fill Dij(m,k) input matrix
                k = mk_mat[x,1]                       # gives the Kth tier of Mth RAT system
                Dij_inp[k-1,0]=mk_mat[x,3]/mk_mat[i,3]
                Dij_inp[k-1,1]=Tmk[x]/Tmk[i]
                Dij_inp[k-1,2]=mk_mat[x,5]
                if(mk_mat[x,2]):
                    Dij_inp[k-1,3]=mk_mat[x,4]
                else:
                    Dij_inp[k-1,4]=mk_mat[x,4]
            Dij_poly = D_ij(Dij_inp,t_x((rate_threshold/bw)*Nij[i]))
            alpha_norm_dij = Dij_inp[:,2]/mk_mat[i,5]
            alpha_norm_dij[alpha_norm_dij==0]=1  #to avoid divide by zero warning :):):)
            f1 = lambda y: sum(Dij_poly*(y**(2/alpha_norm_dij)))
            f2 = lambda y: sum(np.asarray(Gij_mk_array[i])*(y**(2/np.asarray(alpha_norm_mk_array[i]))))
            f3 = lambda y: (t_x((rate_threshold/bw)*Nij[i])*y**(mk_mat[i,5])*0)/mk_mat[i,3]
            f4 = lambda y: -f3(y)-mp.pi*f2(y)-mp.pi*f1(y)
            temp_integral_result[i] = mp.quad(lambda y: y * mp.exp(f4(y)),[0,mp.inf])
        rate_coverage=(sum((2*np.pi*mk_mat[:,4])*temp_integral_result))
    return (rate_coverage)

#-----------------------------------------------------------------------------------------------------------------------