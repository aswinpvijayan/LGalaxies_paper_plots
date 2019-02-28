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
import seaborn as sns
sns.set_context("paper")
sns.set_style(style='white') 

"""
    Plot of depletions as a function of (1) stellar mass (2) metallicity of the corresponding galaxy
"""

redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])

snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])

sims = np.array(['MR', 'MRII'])

inp = int(sys.argv[1])
i = int(sys.argv[2])

filesMR = '../Dust_output/MR/SA_output_5.h5'
filesMRII = '../Dust_output/MRII_sub/SA_output_*'
files = np.array([filesMR, filesMRII])

colors = ['red', 'blue', 'orange', 'green', 'violet', 'indigo', 'violet', 'brown', 'yellow']
arr = np.array(['C', 'O', 'Mg', 'Si', 'Ca', 'Fe'])


def outputs(x, y, sSFR, Type, z):
    
    """
    Function remove_ is for removing any nan's, inifnities and other non-physical values, can specify if the zero values in the array should not be removed, by giving the array number as argument. It then selects the central halos which are above a particular sSFR. Returns the number density in a pixel and the median (in 13 bins in the x axis in logarithmic space) with the values of the 16th percentile and the 84th percentile.
    """
    
    out = get_.remove_(np.array([x, y, sSFR, Type]), np.array([3]))  
    out = out[(out[3] == 0)]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    
    xx, yy, yy_up, yy_low = get_.get_median(thisx, thisy, n = 10)
    
    den = get_.get_density(thisx, thisy)
    
    return thisx, thisy, xx, yy, yy_up, yy_low, den


def plot_Mstar_dep(files, z, axs, snaps, i, on):

    add = sims[i]
    snap = snaps[i][np.where(redshift == str(z))[0][0]]
    Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
    if on and i == 0:
        ok = np.where(Mstar >= 10**9.0)[0]
    elif on and i == 1:
        ok = np.logical_and(Mstar > 10**7.2, Mstar < 10**9.0)
    else:
        ok = np.array([True]*len(Mstar))
    
    Mstar = Mstar[ok]
    Type = get_.get_var(files[i], 'Type', snap)[ok]

    Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
    Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
    Mdust = Mdust1 + Mdust2    
    
    Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok]
    Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]
    Mcg = Mcg1 + Mcg2
    
    SFR = get_.get_var(files[i], 'Sfr', snap)[ok]
    
    sSFR = SFR/Mstar
    ok = np.where(np.nansum(Mcg, axis = 1) > 1e6)
    plot_figure(Mstar[ok], Mcg[ok], Mdust[ok], sSFR[ok], Type[ok], z, axs)
    
    axs.set_xlim((7.5,11.6))
    axs.set_xticks([8, 9, 10, 11])
    
    return add


def plot_Z_dep(files, z, axs, snaps, i, on):

    if on:
        add = sims[0]+'_'+sims[1]
        snap = snaps[0][np.where(redshift == str(z))[0][0]]
        Mstar = (get_.get_var(files[0], 'StellarMass', snap)*1e10)/0.673
        Type = get_.get_var(files[0], 'Type', snap)

        Mdust1 = get_.get_var(files[0], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[0], 'DustColdGasClouds_elements', snap)
        Mdust = Mdust1 + Mdust2    
        
        Mcg1 = get_.get_var(files[0], 'ColdGasDiff_elements', snap)
        Mcg2 = get_.get_var(files[0], 'ColdGasClouds_elements', snap)
        Mcg = Mcg1 + Mcg2
        
        SFR = get_.get_var(files[0], 'Sfr', snap)
        
        
        snap = snaps[1][np.where(redshift == str(z))[0][0]]
        Mstar = np.append(Mstar, (get_.get_var(files[1], 'StellarMass', snap)*1e10)/0.673)
        Type = np.append(Type, get_.get_var(files[1], 'Type', snap))

        tmp1 = get_.get_var(files[1], 'DustColdGasDiff_elements', snap)
        tmp2 = get_.get_var(files[1], 'DustColdGasClouds_elements', snap)
        Mdust = np.append(Mdust, tmp1 + tmp2, axis = 0)
        
        tmp1 = get_.get_var(files[1], 'ColdGasDiff_elements', snap)
        tmp2 = get_.get_var(files[1], 'ColdGasClouds_elements', snap)
        Mcg = np.append(Mcg, tmp1 + tmp2, axis = 0)
        
        SFR = np.append(SFR, get_.get_var(files[1], 'Sfr', snap))
        
    else:
        
        Type = get_.get_var(files[i], 'Type', snap)[ok]

        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
        Mdust = Mdust1 + Mdust2    
        
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok]
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]
        Mcg = Mcg1 + Mcg2
        
        SFR = get_.get_var(files[i], 'Sfr', snap)[ok]
        
    
    sSFR = SFR/Mstar
    
    pltx = np.nansum(Mcg[:,2:], axis = 1)/np.nansum(Mcg, axis = 1)
    
    plot_figure(pltx, Mcg, Mdust, sSFR, Type, z, axs)
    
    axs.set_xlim((-3.6,0.1))
    axs.set_xticks([-3.5, -2.5, -1.5, -0.5])
    
    return add

    
def plot_figure(pltx, Mcg, Mdust, sSFR, Type, z, axs):
    
    for l, k in enumerate([4]):
        
        plty = Mdust[:,k]/Mcg[:,k]
        ok = np.logical_and(np.isnan(plty), np.isinf(plty))
        x, y, xx, yy, yy_up, yy_low, den = outputs(pltx[~ok], plty[~ok], sSFR[~ok], Type[~ok], z)
        x = np.log10(x)
        gridsize = np.array([(int(max(x)-min(x))/0.05), int((max(y)-min(y))/0.05)]).astype(int)    
        axs.hexbin(x, y, gridsize=gridsize, cmap = plt.cm.get_cmap('gist_yarg'), mincnt = 5)
        
        axs.plot(np.log10(xx), yy, color = colors[l], label = arr[l])
        axs.plot(np.log10(xx), yy_up, ls = 'dashed', color = colors[l])
        axs.plot(np.log10(xx), yy_low, ls = 'dashed', color = colors[l])
    
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(18)
        
    axs.set_ylim((-0.05,1.2))
    axs.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0])
    axs.grid(True)
    


fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize=(15, 12), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

if inp == 0:
    
    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(M_{dust}/(M_{\odot}))$'
    savename = 'Depletion_Stellar_'
    
    for z in range(0, 9):
        if i <= 1:
            add = plot_Mstar_dep(files, z, axs[z], snaps, i, False)
        else:
            add1 = plot_Mstar_dep(files, z, axs[z], snaps, 0, True)
            if z == 8:
                axs[z].legend(frameon=False, fontsize=15)
            add2 = plot_Mstar_dep(files, z, axs[z], snaps, 1, True)
            add = add1+'_'+add2
        
            
        axs[z].text(8.2, 1, r'$z = {}$'.format(z), fontsize = 18)
        
elif inp == 1:

    xlab = r'$12+log_{10}(O/H)$'
    ylab = r'$log_{10}(M_{Dust}/(M_{\odot}))$'
    savename = 'Depletion_metallicity_'
    
    for z in range(0, 12):
        if i <= 1:
            add = plot_Z_dep(files, z, axs[z], snaps, i, False)
        else:
            add = plot_Z_dep(files, z, axs[z], snaps, 5, True)
        if z == 11:
            axs[z].legend(frameon=False, fontsize=15)
        axs[z].text(-2.75, 1, r'$z = {}$'.format(z), fontsize = 18)


fig.tight_layout()    
fig.subplots_adjust(bottom=0.09, left = 0.08, wspace=0, hspace=0)
fig.text(0.03, 0.5, r'Fraction', va='center', rotation='vertical', fontsize=22)
fig.text(0.5, 0.04, r'$\mathrm{log}_{10}(\mathrm{M}_{*}/\mathrm{M}_{\odot})$', va='center', fontsize=22)
plt.savefig(savename+add+'_subsample.png')
plt.close()
