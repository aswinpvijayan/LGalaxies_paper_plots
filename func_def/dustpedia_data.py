import numpy as np
from uncertainties import unumpy

themis = 'dustpedia_cigale_results_final_version'
dl14 = 'dustpedia_cigale_results_dl14_final_version'
data = np.genfromtxt('./Obs_data/'+dl14+'.csv', skip_header = 1, delimiter = ',', filling_values = np.nan, usecols = (3,4,17,18))

Mstar, Mstar_err, Mdust, Mdust_err = data[:,0], data[:,1], data[:,2], data[:,3]

Mstar_all = unumpy.uarray(Mstar, Mstar_err)
Mdust_all = unumpy.uarray(Mdust, Mdust_err)

Mstar_all = unumpy.log10(Mstar_all)
Mdust_all = unumpy.log10(Mdust_all)

hubbletype = np.genfromtxt('./Obs_data/hubbletype_dustpedia.txt')
ok = np.logical_and(hubbletype >= -6, hubbletype <= 10)
hubbletype = hubbletype[ok]
Mstar_all = Mstar_all[ok]
Mdust_all = Mdust_all[ok]

