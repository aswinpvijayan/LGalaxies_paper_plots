#!/usr/bin/env python

import sys
import os
import warnings
#if not sys.warnoptions:
#    warnings.simplefilter("ignore")
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
sys.path.append('/lustre/scratch/astro/ap629/my_lib/')
from essent import *

data_Mstar= np.genfromtxt('/lustre/scratch/astro/ap629/Dust_paper_plots/Obs_data/Catinelli_2013_Mstar.txt')
data_MHI = np.genfromtxt('/lustre/scratch/astro/ap629/Dust_paper_plots/Obs_data/Catinelli_2013_MHI.txt')

id_Mstar = data_Mstar[:,0]
id_MHI = data_MHI[:,0]

sel_Mstar = np.where(np.in1d(id_Mstar, id_MHI))[0]
sel_MHI = np.where(np.in1d(id_MHI, id_Mstar))[0]

Mstar = data_Mstar[sel_Mstar][:,4]
MHI = data_MHI[sel_MHI][:,-4]

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(10, 8), sharex=True, sharey=True, facecolor='w', edgecolor='k')

axs.scatter(Mstar, MHI, color = 'red')

file_req = '/lustre/scratch/astro/ap629/test_rmol/rmol_{}/MR/hdf5/lgal_z{}_*.hdf5'.format('clay_f', 0)

Mstar = get_data1('StellarMass', file_req)*(1e10)/h
SFR = get_data1('Sfr', file_req)
sSFR = SFR/Mstar
Type = get_data1('Type', file_req)
ok = np.logical_and(sSFR > sSFR_cut(0), Type == 0)
rmolgal = get_data1('t_exchange_eff', file_req)
Mcg = get_data1('ColdGas_elements', file_req)
rmolgal = get_data1('t_exchange_eff', file_req)
MH = Mcg[:,0]
MHI = MH/(1+rmolgal)


axs.scatter(np.log10(Mstar[ok]), np.log10(MHI[ok]), color = 'blue', alpha = 0.05)
axs.set_xlim((8,12))
plt.show()
