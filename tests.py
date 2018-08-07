#!/usr/bin/env python

import sys
import os
import warnings
#if not sys.warnoptions:
#    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import get_

redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])

snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])

def outputs(x, y, sSFR, Type, z):
    
    ok = (sSFR > get_.sSFR_cut(z)) * (Type == 0)
    thisx = x[ok]
    thisy = y[ok]
    
    return thisx, thisy
    
files = '/lustre/scratch/astro/ap629/test_rmol/dust_pop/MR/SA_output*'  
fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(18, 10), sharex=True, sharey=False, facecolor='w', edgecolor='k')  
axs = axs.ravel()
color = ['red', 'green', 'blue', 'orange', 'violet', 'indigo', 'violet', 'brown', 'yellow']
for j, i in enumerate(range(0,1)):    
    snap = snaps[0][np.where(redshift == str(i))[0][0]]
    Type = get_.get_var(files, 'Type', snap)
    Mstar = (get_.get_var(files, 'StellarMass', snap)*1e10)/0.673

    Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)
    Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)
    Mcg = Mcg1 + Mcg2

    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust = Mdust1 + Mdust2    

    SFR = get_.get_var(files, 'Sfr', snap)
    sSFR = SFR/Mstar
    
    fi = get_.get_var(files, 'f_i', snap)
    fc = get_.get_var(files, 'f_c', snap)
    tacc = get_.get_var(files, 't_acc', snap)
    mu_ = get_.get_var(files, 'mu_', snap)
    mu_gas = get_.get_var(files, 'mu_gas', snap)
    
    x, y = outputs(np.log10(Mstar), np.log10(fi), sSFR, Type, i)
    axs[0].scatter(x, y, marker = 'o', s=2)
    axs[0].grid()
    axs[0].set_xlim((7,12))
    axs[0].set_xlabel(r'$log_{10}(M_*)$', fontsize = 16)
    axs[0].set_ylabel(r'$log_{10}(fi)$', fontsize = 16)
    
    x, y = outputs(np.log10(Mstar), np.log10(fc), sSFR, Type, i)
    axs[1].scatter(x, y, marker = 'o', s=2)
    axs[1].grid()
    axs[1].set_xlim((7,12))
    axs[1].set_xlabel(r'$log_{10}(M_*)$', fontsize = 16)
    axs[1].set_ylabel(r'$log_{10}(fc)$', fontsize = 16)
    
    x, y = outputs(np.log10(Mstar), np.log10(tacc), sSFR, Type, i)
    axs[2].scatter(x, y, marker = 'o', s=2)
    axs[2].grid()
    axs[2].set_xlim((7,12))
    axs[2].set_xlabel(r'$log_{10}(M_*)$', fontsize = 16)
    axs[2].set_ylabel(r'$log_{10}(t_{acc})$', fontsize = 16)
    
    x, y = outputs(np.log10(Mstar), np.log10(mu_), sSFR, Type, i)
    axs[3].scatter(x, y, marker = 'o', s=2)
    axs[3].grid()
    axs[3].set_xlim((7,12))
    axs[3].set_xlabel(r'$log_{10}(M_*)$', fontsize = 16)
    axs[3].set_ylabel(r'$log_{10}(\mu_{cloud})$', fontsize = 16)
    
    x, y = outputs(np.log10(Mstar), np.log10(mu_gas), sSFR, Type, i)
    axs[4].scatter(x, y, marker = 'o', s=2)
    axs[4].grid()
    axs[4].set_xlim((7,12))
    axs[4].set_xlabel(r'$log_{10}(M_*)$', fontsize = 16)
    axs[4].set_ylabel(r'$log_{10}(\mu_{gas})$', fontsize = 16)
    
    x, y = outputs(np.log10(Mstar), np.log10(sSFR), sSFR, Type, i)
    axs[5].scatter(x, y, marker = 'o', s=2)
    axs[5].grid()
    axs[5].set_xlim((7,12))
    axs[5].set_xlabel(r'$log_{10}(M_*)$', fontsize = 16)
    axs[5].set_ylabel(r'$log_{10}(sSFR)$', fontsize = 16)

fig.suptitle('GK11', fontsize=16)
#axs.legend(frameon=False)
#axs.set_xlim((8,12))
#axs.grid()
plt.show()
