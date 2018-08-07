#!/usr/bin/env python

import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.append('/lustre/scratch/astro/ap629/my_lib/')
from essent import *
from joblib import Parallel, delayed
import timeit


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):

    return [atoi(c) for c in re.split('(\d+)', text)]

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):

    return [atoi(c) for c in re.split('(\d+)', text)]

def _files(file_req):
    
    if np.isscalar(file_req):
        files = glob.glob(file_req)
        files.sort(key=natural_keys)
    else:
        files = []
        for i in file_req:
            tmp = glob.glob(i)
            tmp.sort(key=natural_keys)
            files.append(tmp)
        files = np.concatenate(files)
    return files

def find(obj, file_req, snap):
    
    start = timeit.default_timer()
    files = _files(file_req)
    files = files[98:100]
    print ('Reading {} from snap {}'.format(obj, snap))
    results = 0.
    for filename in files:
        g = h5py.File(filename, 'r')
        f_0 = g[snap]
        if np.isscalar(results):
            #print ('Reading in {} from {}...'.format(obj, files[0]))
            out = np.array(f_0[obj])
            tmp = np.shape(out)
            if len(tmp) > 1:
                results = np.zeros((int(10**8), tmp[1]))*np.nan
            else:
                results = np.zeros(int(10**8))*np.nan
            length = len(out)
            results[:length] = out
        else:
            out = np.array(f_0[obj])
            tmp = len(out) + length
            results[length:tmp] = out
            length = tmp
        g.close()   
    
    stop = timeit.default_timer()
    print ('Time taken is {}s'.format(stop - start)) 
    
    return out


def obtain(obj, file_req, snap):
    
    start = timeit.default_timer()
    
    #print ('Reading in {} from snap {}...'.format(obj, snap))
    #filess = get_files(files)
    #print (filess)
    out = Parallel(n_jobs = 4, backend="multiprocessing", batch_size = 1)(delayed(find)(i, file_req, snap) for i in obj)   
  
    stop = timeit.default_timer()
    print ('Time taken is {}s'.format(stop - start)) 
        
    return out






def outputs(x, y, sSFR, Type, z):
    
    """
    #Function remove_ is for removing any nan's, inifnities and other non-physical values, can specify if the zero values in the array should not be removed, by giving the array number as argument. It then selects the central halos which are above a particular sSFR. Returns the number density in a pixel and the median (in 13 bins in the x axis in logarithmic space) with the values of the 16th percentile and the 84th percentile.
    """
    
    out = remove_(np.array([x, y, sSFR, Type]), np.array([3]))  
    out = out[(out[2] > sSFR_cut(z)) & (out[3] == 0)]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    this_sSFR = out[2]
    this_Type = out[3]
    
    del out, x, y, sSFR, Type
    
    xx, yy, yy_up, yy_low = get_median(thisx, thisy, n = 25)
    
    den = get_density(thisx, thisy)
    
    return thisx, thisy, den, xx, yy, yy_up, yy_low


fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize=(15, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

redshift = np.array(['0', '0.11', '0.26', '0.52', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
snapnum = ['58', '54', '50', '45', '38', '30', '25', '22', '19', '17', '15', '13', '12']

xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
ylab = r'$log_{10}(M_{Dust}/(M_{\odot}))$'
savename = 'Dust_Stellar_MR'

#files = '/lustre/scratch/astro/ap629/test_rmol/dust_new/MR/SA_output_5.h5'
files = '/lustre/scratch/astro/ap629/Dust_output_22May/MR/SA_output_*'
#f = h5py.File(files, 'r')
vars_ = np.array(['Type', 'StellarMass', 'DustColdGasDiff_elements', 'DustColdGasClouds_elements', 'Sfr'])
from obs_plots import DM_obs
for z in range(0, 9):
    
    snap = snapnum[np.where(redshift == str(z))[0][0]]
    #g = f[snap]
    results = obtain(vars_, files, snap)
    Type, Mstar, DustColdGasDiff_elements, DustColdGasClouds_elements, Sfr = results
    Mstar = (Mstar*1e10)/0.673
    Mdust1 = np.sum(DustColdGasDiff_elements, axis = 1)  
    Mdust2 = np.sum(DustColdGasClouds_elements, axis = 1) 
    Mdust = Mdust1 + Mdust2
    sSFR = Sfr/Mstar
    
    x, y, den, xx, yy, yy_up, yy_low = outputs(Mstar, Mdust, sSFR, Type, z)
    
    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    xlim = [6.5,11.5]
    ylim = [2,9.95]
    xticks = [7, 8, 9, 10, 11, 12]
    axs[z].text(8.75, 9, r'$z = {}$'.format(z), fontsize = 18)
    
    DM_obs(axs[z], z)   #Plotting the observational data points

    axs[z].scatter(x, y, marker = 'o', s=10, c = den, edgecolors='None', alpha = 0.08, cmap = plt.cm.get_cmap('gist_gray'))
    
    axs[z].plot(xx, yy, lw = 2, color = 'grey')
    axs[z].plot(xx, yy_up, lw = 2, ls = 'dashed', color = 'grey')
    axs[z].plot(xx, yy_low, lw = 2, ls = 'dashed', color = 'grey')
    
    del x, y, xx, yy, yy_up, yy_low
    
    axs[z].set_xticks(xticks)
    
    for label in (axs[z].get_xticklabels() + axs[z].get_yticklabels()):
        label.set_fontsize(18)
    
    axs[z].grid()
    lgd = axs[z].legend(frameon=False, fontsize = 16, markerscale=2, loc = 4, numpoints=1, handletextpad=0.005)
    if np.isscalar(lgd):
        lgd.set_zorder(100)
    axs[z].set_xlim(xlim)
    axs[z].set_ylim(ylim)
    
fig.tight_layout()    
fig.subplots_adjust(bottom=0.09, left = 0.08, wspace=0, hspace=0)
fig.text(0.03, 0.5, ylab, va='center', rotation='vertical', fontsize=22)
fig.text(0.5, 0.04, xlab, va='center', fontsize=22)
plt.savefig('Mstar_Mdust1.png')
plt.show()

