#!/usr/bin/env python
import datetime
print ('Time at the start is: ', str(datetime.datetime.now()))

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
import mbb 
import gc
import seaborn as sns
sns.set_context("paper")
from astropy.cosmology import Planck13

        ###################################################################################
h = 0.673        
redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])
sims = np.array(['MR', 'MRII'])

MR_vol = (480.279/h)**3  #Millennium
MRII_vol = (96.0558/h)**3 #Millennium II
fac = MR_vol/MRII_vol

colors = ['orange', 'violet', 'cyan', 'magenta', 'brown', 'orange', 'black']

def outputs(x, y, sSFR, Type, z):
    
    """
    #Function remove_ is for removing any nan's, inifnities and other non-physical values, can specify if the zero values in the array should not be removed, by giving the array number as argument. It then selects the central halos which are above a particular sSFR. Returns the number density in a pixel and the median with the values of the 16th percentile and the 84th percentile.
    """
    
    out = get_.remove_(np.array([x, y, sSFR, Type]), np.array([3]))  
    out = out[(out[2] > get_.sSFR_cut(z)) & (out[3] == 0)]
    #out = out[out[3] == 0]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    thissSFR = out[2]
    
    del out, x, y, sSFR, Type
    
    xx, yy, yy_up, yy_low = get_.get_median(thisx, thisy, n = 10)
    
    return xx, yy, yy_up, yy_low


def plot_Mstar_Mdust(fil, z, axs, snapnum, i, on, k, label):
    
    add = sims[i] 
    snap = snapnum[i][np.where(redshift == str(z))[0][0]]
    print (fil[i])
    Mstar = (get_.get_var(fil[i], 'StellarMass', snap)*1e10)/0.673
    if on and i == 0:
        ok = np.where(Mstar >= 10**8.9)[0]
    elif on and i == 1:
        ok = np.logical_and(Mstar > 10**7.5, Mstar < 10**8.9)
    else:
        ok = np.array([True]*len(Mstar))
        
    Mstar = Mstar[ok]
    Type = get_.get_var(fil[i], 'Type', snap)[ok]
    Mdust1 = get_.get_var(fil[i], 'DustColdGasDiff_elements', snap)[ok]
    Mdust2 = get_.get_var(fil[i], 'DustColdGasClouds_elements', snap)[ok]
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    Mdust = Mdust1 + Mdust2
    SFR = get_.get_var(fil[i], 'Sfr', snap)[ok]
    sSFR = SFR/Mstar
    
    xx, yy, yy_up, yy_low = outputs(Mstar, Mdust, sSFR, Type, z)
    
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    plot_figure(axs, z, xx, yy, yy_up, yy_low, k, label)
    
    xlim = [7.5,11.7]
    ylim = [1.5,10.5]
    xticks = [8, 9, 10, 11]
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    
    axs.set_xticks(xticks)
    
    return add


def plot_figure(axs, z, xx, yy, yy_up, yy_low, k, label):
    
    if label != None:
        axs.plot(xx, yy, lw = 1, color = colors[k], label = r'$t_{acc}=$'+label+r'$yrs$')
    else:
        axs.plot(xx, yy, lw = 1, color = colors[k])
    axs.plot(xx, yy_up, lw = 1, ls = 'dashed', color = colors[k])
    axs.plot(xx, yy_low, lw = 1, ls = 'dashed', color = colors[k])
    
    del xx, yy, yy_up, yy_low
    
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(18)
    
    axs.grid(True)
    
i = int(sys.argv[1])        

filesMR = np.array(['../test_rmol/def_1e3/MR/SA_output_*', '../test_rmol/def_5e3/MR/SA_output_*', '../ap629/test_rmol/def_15e3/MR/SA_output_*', '../test_rmol/def_1e5/MR/SA_output_*', '../test_rmol/def_1e6/MR/SA_output_*'])
filesMRII = np.array(['../test_rmol/def_1e3/MRII/SA_output_*', '../test_rmol/def_5e3/MRII/SA_output_*', '../test_rmol/def_15e3/MRII/SA_output_*', '../test_rmol/def_1e5/MRII/SA_output_*', '../test_rmol/def_1e6/MRII/SA_output_*'])


files = np.hstack([np.array([filesMR]).T, np.array([filesMRII]).T])

labels = np.array([r'$10^3$', r'$5\times 10^3$', r'$15\times 10^3$', r'$10^5$', r'$10^6$'])
titles = np.array(['1e3', '5e3', '15e3', '1e5', '1e6'])

fig, axs = plt.subplots(nrows = 1, ncols = 4, figsize=(20, 8), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
ylab = r'$log_{10}(M_{dust}/(M_{\odot}))$'
savename = 'Dust_Stellar_diff_taccs'

from obs_plots import DM_obs 
z_ = np.array([0, 2, 4, 6])
for j, z in enumerate(z_):
    DM_obs(axs[j], z)   #Plotting the observational data points
    for k in range(0, len(labels)): 
        if i <= 1:
            add = plot_Mstar_Mdust(files[k], z, axs[j], snaps, i, False, k, labels[k])
            lgd = axs[j].legend(fontsize = 16, markerscale=2, loc = 4, numpoints=1, handletextpad=0.005)
            if np.isscalar(lgd):
                lgd.set_zorder(100)
        else:
            if z == 6:
                add1 = plot_Mstar_Mdust(files[k], z, axs[j], snaps, 0, True, k, labels[k])
            else:
                add1 = plot_Mstar_Mdust(files[k], z, axs[j], snaps, 0, True, k, None)
            lgd = axs[j].legend(fontsize = 16, markerscale=2, loc = 4, numpoints=1, handletextpad=0.005, frameon = False)
            if np.isscalar(lgd):
                lgd.set_zorder(100)
            add2 = plot_Mstar_Mdust(files[k], z, axs[j], snaps, 1, True, k, None)
            add = add1+'_'+add2
            
    axs[j].text(8.65, 9.5, r'$z = {}$'.format(z), fontsize = 18)
    
    
fig.tight_layout()    
fig.subplots_adjust(bottom=0.1, left = 0.08, wspace=0, hspace=0)
fig.text(0.03, 0.5, ylab, va='center', rotation='vertical', fontsize=22)
fig.text(0.5, 0.03, xlab, va='center', fontsize=22)
plt.savefig(savename+add+'_taccs.eps')
plt.close()
print ('End time is ', str(datetime.datetime.now()))
