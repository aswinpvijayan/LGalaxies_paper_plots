import numpy as np
from uncertainties import unumpy
  
data = np.genfromtxt('./Obs_data/De_Cia_2016.txt', delimiter=',')
name, DTM, DTMerr, z, OH, OHerr, xuplims, yuplims = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7]

DTM_all = unumpy.log10(unumpy.uarray(DTM, DTMerr)*0.464)
DTM_all = 0.98*DTM_all/(0.98*DTM_all+1.)

OH_abund = unumpy.uarray(OH, OHerr) + 8.69

z = np.round(z,0)


def delta_O(ZnFe):
    
    A = unumpy.uarray([-0.02], [0.1])
    B = unumpy.uarray([-0.15], [0.09])
    
    return A + B*ZnFe

data = np.genfromtxt('./Obs_data/De_Cia_a_2016.txt', delimiter=',')
name, ZnFe, ZnFe_err, MH, MH_err, DTM, DTMerr, yuplims_a, z_a, xuplims_a = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8], data[:,9]

tmp = unumpy.log10(unumpy.uarray(DTM, DTMerr))
DTM_all = np.append(DTM_all, unumpy.log10(unumpy.uarray(DTM, DTMerr)*0.464))

del_O = delta_O(unumpy.uarray(ZnFe, ZnFe_err))
OH_abund = np.append(OH_abund, del_O + unumpy.uarray(MH, MH_err) + 8.69)

z = np.append(z, np.round(z_a,0))
xuplims = np.append(xuplims, xuplims_a)
yuplims = np.append(yuplims, yuplims_a)

