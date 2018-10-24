#!/usr/bin/env python

import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
sys.path.append('func_def/')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import get_

redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])

snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snapnum = np.array([snapnumMR, snapnumMRII])

filesMR = '../Rob_dust_output/MR/SA_output_*'
filesMRII = '../Rob_dust_output/MRII/SA_output_*'
files = np.array([filesMR, filesMRII])
fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize=(10, 8), sharex=False, sharey=True, facecolor='w', edgecolor='k')  
axs = axs.ravel()
color = ['red', 'green', 'blue', 'orange', 'violet', 'indigo', 'violet', 'brown', 'yellow']
for j, k in enumerate(range(6,7)):    
    
    i = 0
    z = 6
    snap = snapnum[i][np.where(redshift == str(z))[0][0]]
    Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
    ok = np.where(Mstar >= 10**8.85)[0]
    Mstar = Mstar[ok]
    Type = get_.get_var(files[i], 'Type', snap)[ok]
    dustrates = get_.get_var(files[i], 'DustColdGasRates', snap)[ok]
    Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok]
    Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]
    Mcg = Mcg1 + Mcg2
    
    Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
    Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
    Mdust = Mdust1 + Mdust2
    
    met = (Mcg[:,4]-Mdust[:,4])/(15.9994*(Mcg[:,0]))  #O/H ratio
    
    Mcg = np.sum(Mcg[:,2:], axis = 1) #Total cold gas mass
    
    Mdust = np.nansum(Mdust, axis = 1) 
    
    Mratio = Mdust/Mcg   #Dust-to-total metal mass ratio. Dimensionless.
    
    i = 1
    snap = snapnum[i][np.where(redshift == str(z))[0][0]]
    tmp = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
    ok = np.logical_and(tmp > 10**5.95, tmp < 10**8.9)
    Mstar = np.append(Mstar, tmp[ok])
    Type = np.append(Type, get_.get_var(files[i], 'Type', snap)[ok])
    dustrates = np.append(dustrates, get_.get_var(files[i], 'DustColdGasRates', snap)[ok], axis = 0)
    
    Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok]
    Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]
    tmp1 = Mcg1 + Mcg2
    
    Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
    Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
    tmp2 = Mdust1 + Mdust2
    
    met = np.append(met, (tmp1[:,4]-tmp2[:,4])/(15.9994*tmp1[:,0]))
    
    tmp1 = np.nansum(tmp1[:,2:], axis = 1)
    Mcg = np.append(Mcg, tmp1)
    
    tmp2 = np.nansum(tmp2, axis = 1)
    Mdust = np.append(Mdust, tmp2) 
    
    Mratio = np.append(Mratio, tmp2/tmp1)
    
ok = np.logical_and(Mdust > 0.0, np.logical_and(12. + np.log10(met) > 6.4, Type == 0))

Mstar = Mstar[ok]
Mratio = Mratio[ok]
Mdust = Mdust[ok]
met = met[ok]
dustrates = dustrates[ok]

ok = np.logical_and(np.log10(Mratio) > -2, np.log10(Mratio) > -2)

xx, yy, yy_up, yy_low = get_.get_median(Mstar[ok], dustrates[:,1][ok], n = 16)
axs[0].plot(np.log10(xx), np.log10(yy), color = 'red', label = r'$SNII$')
axs[0].plot(np.log10(xx), np.log10(yy_up), lw = 1, ls = 'dashed', color = 'red')
axs[0].plot(np.log10(xx), np.log10(yy_low), lw = 1, ls = 'dashed', color = 'red')
xx, yy, yy_up, yy_low = get_.get_median(Mstar[ok], dustrates[:,3][ok], n = 16)
axs[0].plot(np.log10(xx), np.log10(yy), color = 'green', label = r'$GG$')
axs[0].plot(np.log10(xx), np.log10(yy_up), lw = 1, ls = 'dashed', color = 'green')
axs[0].plot(np.log10(xx), np.log10(yy_low), lw = 1, ls = 'dashed', color = 'green')
axs[0].set_xlim((7.5, 10.7))
axs[0].set_xticks((8, 9, 10))
axs[0].set_ylim((-7, 1.5))
axs[0].grid(True)

xx, yy, yy_up, yy_low = get_.get_median(met[ok], dustrates[:,1][ok], n = 13)
axs[1].plot(12. + np.log10(xx), np.log10(yy), color = 'red', label = r'$SNII$')
axs[1].plot(12. + np.log10(xx), np.log10(yy_up), lw = 1, ls = 'dashed', color = 'red')
axs[1].plot(12. + np.log10(xx), np.log10(yy_low), lw = 1, ls = 'dashed', color = 'red')
xx, yy, yy_up, yy_low = get_.get_median(met[ok], dustrates[:,3][ok], n = 13)
axs[1].plot(12. + np.log10(xx), np.log10(yy), color = 'green', label = r'$GG$')
axs[1].plot(12. + np.log10(xx), np.log10(yy_up), lw = 1, ls = 'dashed', color = 'green')
axs[1].plot(12. + np.log10(xx), np.log10(yy_low), lw = 1, ls = 'dashed', color = 'green')
axs[1].set_xlim(([6.5,9.7]))
axs[1].set_xticks((7, 8, 9))
axs[1].set_ylim((-7, 1.5))
axs[1].grid(True)

xx, yy, yy_up, yy_low = get_.get_median(Mstar[~ok], dustrates[:,1][~ok], n = 16)
axs[2].plot(np.log10(xx), np.log10(yy), color = 'red', label = r'$SNII$')
axs[2].plot(np.log10(xx), np.log10(yy_up), lw = 1, ls = 'dashed', color = 'red')
axs[2].plot(np.log10(xx), np.log10(yy_low), lw = 1, ls = 'dashed', color = 'red')
xx, yy, yy_up, yy_low = get_.get_median(Mstar[~ok], dustrates[:,3][~ok], n = 16)
axs[2].plot(np.log10(xx), np.log10(yy), color = 'green', label = r'$GG$')
axs[2].plot(np.log10(xx), np.log10(yy_up), lw = 1, ls = 'dashed', color = 'green')
axs[2].plot(np.log10(xx), np.log10(yy_low), lw = 1, ls = 'dashed', color = 'green')
axs[2].set_xlim((7.5, 10.7))
axs[2].set_xticks((8, 9, 10))
axs[2].set_ylim((-7, 1.5))
axs[2].grid(True)
axs[2].set_xlabel(r'$log_{10}(M_{*}/M_{\odot})$', fontsize = 22)

xx, yy, yy_up, yy_low = get_.get_median(met[~ok], dustrates[:,1][~ok], n = 16)
axs[3].plot(12. + np.log10(xx), np.log10(yy), color = 'red', label = r'$SNII$')
axs[3].plot(12. + np.log10(xx), np.log10(yy_up), lw = 1, ls = 'dashed', color = 'red')
axs[3].plot(12. + np.log10(xx), np.log10(yy_low), lw = 1, ls = 'dashed', color = 'red')
xx, yy, yy_up, yy_low = get_.get_median(met[~ok], dustrates[:,3][~ok], n = 16)
axs[3].plot(12. + np.log10(xx), np.log10(yy), color = 'green', label = r'$GG$')
axs[3].plot(12. + np.log10(xx), np.log10(yy_up), lw = 1, ls = 'dashed', color = 'green')
axs[3].plot(12. + np.log10(xx), np.log10(yy_low), lw = 1, ls = 'dashed', color = 'green')
axs[3].set_xlim((6.5,9.7))
axs[3].set_ylim((-7, 1.5))
axs[3].set_xticks((7, 8, 9))
axs[3].grid(True)
axs[3].set_xlabel(r'$12+log_{10}(O/H)$', fontsize = 22)

for p in range(0, 4):
    for label in (axs[p].get_xticklabels() + axs[p].get_yticklabels()):
        label.set_fontsize(17)

fig.subplots_adjust(hspace=0, wspace = 0)
fig.text(0.05, 0.5, r'$log_{10}(\Phi_{\mathrm{DR}}/(M_{\odot}\mathrm{yr}^{-1}))$', va='center', rotation='vertical', fontsize=22)
fig.savefig('z6_DTM_low_high.eps')
plt.show()
