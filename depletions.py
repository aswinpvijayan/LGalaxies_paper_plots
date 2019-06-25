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
    Figure 2 in paper
"""

redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])

snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])

sims = np.array(['MR', 'MRII'])

inp = int(sys.argv[1])
i = int(sys.argv[2])

filesMR = '../Dust_output/MR/SA_output_*'
filesMRII = '../Dust_output/MRII/SA_output_*'
files = np.array([filesMR, filesMRII])

colors = ['blue', 'brown', 'orange', 'green', 'violet', 'indigo', 'violet', 'red', 'yellow']
arr = np.array(['C', 'O', 'Mg', 'Si', 'Ca', 'Fe'])
clmaps = ['gist_yarg', 'Reds']

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
        ok = np.where(Mstar>=1e9)
        Mstar = Mstar[ok]
        Type = get_.get_var(files[0], 'Type', snap)[ok]

        Mdust1 = get_.get_var(files[0], 'DustColdGasDiff_elements', snap)[ok]
        Mdust2 = get_.get_var(files[0], 'DustColdGasClouds_elements', snap)[ok]
        Mdust = Mdust1 + Mdust2    
        
        Mcg1 = get_.get_var(files[0], 'ColdGasDiff_elements', snap)[ok]
        Mcg2 = get_.get_var(files[0], 'ColdGasClouds_elements', snap)[ok]
        Mcg = Mcg1 + Mcg2
        
        SFR = get_.get_var(files[0], 'Sfr', snap)[ok]
        
        
        snap = snaps[1][np.where(redshift == str(z))[0][0]]
        tmp = (get_.get_var(files[1], 'StellarMass', snap)*1e10)/0.673
        ok = np.logical_and(tmp > 1e7, tmp <= 1e9)
        Mstar = np.append(Mstar, tmp[ok])
        
        Type = np.append(Type, get_.get_var(files[1], 'Type', snap)[ok])

        tmp1 = get_.get_var(files[1], 'DustColdGasDiff_elements', snap)[ok]
        tmp2 = get_.get_var(files[1], 'DustColdGasClouds_elements', snap)[ok]
        Mdust = np.append(Mdust, tmp1 + tmp2, axis = 0)
        
        tmp1 = get_.get_var(files[1], 'ColdGasDiff_elements', snap)[ok]
        tmp2 = get_.get_var(files[1], 'ColdGasClouds_elements', snap)[ok]
        Mcg = np.append(Mcg, tmp1 + tmp2, axis = 0)
        
        SFR = np.append(SFR, get_.get_var(files[1], 'Sfr', snap)[ok])
        
    else:
        
        add = sims[i]
        snap = snaps[0][np.where(redshift == str(z))[0][0]]
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        if i == 0: 
            ok = np.where(Mstar>=1e9)[0]
        else:
            ok = np.logical_and(Mstar > 1e7, Mstar <= 1e9)
            
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
    
    pltx = (Mcg[:,4] - Mdust[:,4])/(15.9994*(Mcg[:,0]))
    
    plot_figure(pltx, Mcg, Mdust, sSFR, Type, z, axs)
    
    xlim = [6.8,9.8]
    xticks = [7, 7.5, 8, 8.5, 9, 9.5]
    
    axs.set_xlim(xlim)
    axs.set_xticks(xticks)
    
    return add

    
def plot_figure(pltx, Mcg, Mdust, sSFR, Type, z, axs):
    
    for l, k in enumerate([2, 4]):
        
        plty = Mdust[:,k]/Mcg[:,k]
        ok = np.logical_and(np.isnan(plty), np.isinf(plty))
        x, y, xx, yy, yy_up, yy_low, den = outputs(pltx[~ok], plty[~ok], sSFR[~ok], Type[~ok], z)
        x = 12. + np.log10(x)
        gridsize = np.array([(int(max(x)-min(x))/0.03), int((max(y)-min(y))/0.03)]).astype(int)    
        axs.hexbin(x, y, gridsize=gridsize, bins = 'log', cmap = plt.cm.get_cmap(clmaps[l]), mincnt = 1)
        
        axs.plot(12. + np.log10(xx), yy, color = colors[l], label = arr[l])
        axs.plot(12. + np.log10(xx), yy_up, ls = 'dashed', color = colors[l])
        axs.plot(12. + np.log10(xx), yy_low, ls = 'dashed', color = colors[l])
    
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(18)
        
    axs.set_ylim((0.0,0.76))
    axs.set_yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
    axs.legend(frameon=False, fontsize=17, loc='center right')
    axs.grid(True)
    


fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(7, 7), sharex=True, sharey=True, facecolor='w', edgecolor='k')
#axs = axs.ravel()

if inp == 0:
    
    xlab = r'$log_{10}(\mathrm{M}_{*}/(\mathrm{M}_{\odot}))$'
    savename = 'Depletion_Stellar_'
    
    for z in range(0, 9):
        if i <= 1:
            add = plot_Mstar_dep(files, z, axs[z], snaps, i, False)
            axs[z].legend(frameon=False, fontsize = 18)
        else:
            add1 = plot_Mstar_dep(files, z, axs[z], snaps, 0, True)
            axs[z].legend(frameon=False, fontsize = 18)
            if z == 8:
                axs[z].legend(frameon=False, fontsize=15)
            add2 = plot_Mstar_dep(files, z, axs[z], snaps, 1, True)
            add = add1+'_'+add2
        
            
        axs[z].text(8.2, 0.85, r'$z = {}$'.format(z), fontsize = 18)
        
        
elif inp == 1:

    xlab = r'$12 + \mathrm{log}_{10}(\mathrm{O/H})$'
    savename = 'Depletion_metallicity_'
    
    for z in range(0, 1):
        if i <= 1:
            add = plot_Z_dep(files, z, axs, snaps, i, False)
        else:
            add = plot_Z_dep(files, z, axs, snaps, 5, True)
        #axs.text(-1., 0.45, r'$z = {}$'.format(z), fontsize = 18)


fig.tight_layout()    
fig.subplots_adjust(bottom=0.12, left = 0.14, wspace=0, hspace=0)
fig.text(0.015, 0.52, r'$\mathrm{Element\ Depletion\ Fraction}$', va='center', rotation='vertical', fontsize=22)
fig.text(0.45, 0.04, xlab, va='center', fontsize=22)
plt.savefig(savename+add+'_z0.pdf')
plt.close()
