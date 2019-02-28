import numpy as np
from uncertainties import unumpy

def X_H(F, M_H):
    
    #A = unumpy.uarray([-0.02], [0.1])
    #B = unumpy.uarray([-0.15], [0.09])
    
    A = unumpy.uarray([-0.145], [0.051])
    B = unumpy.uarray([-0.225], [0.053])
    
    return A + B*(F-0.598) + M_H
    #return A + B*(F-0.598)
    
data = np.genfromtxt('./Obs_data/Wiseman_data_2017.txt', delimiter=',')
F, Ferr, M_H, M_Herror, z, DTM, DTMerr = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]

F_all = unumpy.uarray(F, Ferr)
M_H_all = unumpy.uarray(M_H, M_Herror)
DTM_all = unumpy.log10(unumpy.uarray(DTM, DTMerr)*0.464)

depl = X_H(F_all, M_H_all)  
OH_abund = depl + 8.69

z = np.round(z,0)
