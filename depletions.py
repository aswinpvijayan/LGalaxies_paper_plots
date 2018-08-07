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
import seaborn as sns
sns.set_context("paper")

redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])

snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])

def outputs(x, y, sSFR, Type, z):
    
    """
    #Function remove_ is for removing any nan's, inifnities and other non-physical values, can specify if the zero values in the array should not be removed, by giving the array number as argument. It then selects the central halos which are above a particular sSFR. Returns the number density in a pixel and the median (in 13 bins in the x axis in logarithmic space) with the values of the 16th percentile and the 84th percentile.
    """
    #out = np.array([x, y, sSFR, Type])
    #out = get_.remove_(np.array([x, y, sSFR, Type]), np.array([3]))  
    ok = (sSFR > get_.sSFR_cut(z)) * (Type == 0)
    thisx = x[ok]
    thisy = y[ok]
    
    #del out, x, y, sSFR, Type
    
    
    #xx, yy, yy_up, yy_low = get_.get_median(thisx, thisy, n = 25)
    
    #den = get_.get_density(thisx, thisy)
    
    #return thisx, thisy, den, xx, yy, yy_up, yy_low
    return thisx, thisy
    
files = '/lustre/scratch/astro/ap629/test_rmol/dust_krumholtz/MR/SA_output*'  
#fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(15, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')  
color = ['red', 'green', 'blue', 'orange', 'violet', 'indigo', 'violet', 'brown', 'yellow']
fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(15, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()
arr = np.array(['C', 'O', 'Mg', 'Si', 'Ca', 'Fe'])
for j, i in enumerate(range(0,1)):    
    snap = snaps[0][np.where(redshift == str(i))[0][0]]
    Type = get_.get_var(files, 'Type', snap)
    mu_gas = get_.get_var(files, 'mu_gas', snap)
    Mstar = (get_.get_var(files, 'StellarMass', snap)*1e10)/0.673

    Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)
    Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)
    Mcg = Mcg1 + Mcg2

    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust = Mdust1 + Mdust2    

    SFR = get_.get_var(files, 'Sfr', snap)
    sSFR = SFR/Mstar
    DGHratio = np.sum(Mdust, axis = 1)/Mcg[:,0]
    
    O_dep = np.log10((Mcg[:,4]-Mdust[:,4])/(15.999*Mcg[:,0])) - np.log10(477./909964.)
    Mg_dep = np.log10((Mcg[:,6]-Mdust[:,6])/(24.305*Mcg[:,0])) - np.log10(36./909964.)
    Si_dep = np.log10((Mcg[:,7]-Mdust[:,7])/(28.085*Mcg[:,0]))- np.log10(30./909964.)
    Ca_dep = np.log10((Mcg[:,9]-Mdust[:,9])/(40.078*Mcg[:,0])) - np.log10(2./909964.)
    Fe_dep = np.log10((Mcg[:,10]-Mdust[:,10])/(55.845*Mcg[:,0])) - np.log10(30./909964.)
    
    Fe_met = np.log10((Mcg[:,10])/(55.845*Mcg[:,0])) - np.log10(30./909964.)
    
    for l, k in enumerate([2,4,6,7,9,10]):
        x, y = outputs(Mstar, (Mdust[:,k])/(Mcg[:,k]), sSFR, Type, i)
        
        axs[l].scatter(np.log10(x), y, marker = 'o', s=5, c = color[l], label = arr[l])
        axs[l].set_xlim((7,11.6))
        axs[l].set_ylim((0,1.1))
        axs[l].grid()
        axs[l].legend(frameon=False, fontsize=15)

fig.subplots_adjust(bottom=0.09, left = 0.08, wspace=0, hspace=0)
#axs.legend(frameon=False)
#axs.set_xlabel(r'$12+log_{10}(O/H)$', fontsize = 25)
axs[-2].set_xlabel(r'$log_{10}(M_{*}/M_{\odot})$', fontsize = 22)
#axs.set_ylabel(r'$\frac{O(dust)}{O(total)}$', fontsize = 22)
#axs.set_ylim((-6,-1))
#axs.set_xlim((6,11))
#axs.grid()
plt.savefig('depletions.png')
plt.show()
