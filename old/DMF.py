#!/usr/bin/env python

import sys
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
sys.path.append('/lustre/scratch/astro/ap629/my_lib/')
from essent import *


def files_(turn_on, z):
    
    """
    Function to pick the files required for the analysis. MRII is added or not, depending on the value of 'turn_on'. 
    """
    
    if turn_on == 1:
        file_req = ['/lustre/scratch/astro/ap629/Dust_output_17jan/MR/hdf5/lgal_z{}_N*.hdf5'.format(z), '/lustre/scratch/astro/ap629/Dust_output_17jan/MRII/hdf5/lgal_z{}_N*.hdf5'.format(z)] 
    else:
        file_req = '/lustre/scratch/astro/ap629/Dust_output_17jan/MR/hdf5/lgal_z{}_N*.hdf5'.format(z)
    
    return file_req

def dust_massfn(file_req, vol):

    """
    Outputs the dust mass function for input files of a particular redshift, for a given volume
    """

    Mdust = get_data1('Dust_elements', file_req)
    Mdust = np.sum(Mdust, axis = 1)
    out = remove_(Mdust)
    out = np.array(out[0])
    #out = out[out > 1e5]
    #print (min(out), max(out))
    dbins = np.arange(np.log10(min(out))-0.2, np.log10(max(out))+0.2, 0.2)
    bins, edges = np.histogram(np.log10(out), dbins)
    
    xx = (edges[1:]+edges[:-1])/2
    binsize = (xx[-1] - xx[0])/len(xx)
    yy = bins/(vol*binsize)   
    
    return xx, yy 


h = 0.673 #little h as used in the model

#Simulation volumes
MR_vol = (480.279/h)**3  #Millennium
MRII_vol = (96.0558/h)**3 #Millennium II

#turn_on = 0 #Selecting MR and MRII, turn_on = 1 selects MRII for plotting as well.
turn_on = int(input('Enter 1 for plotting the dust mass function of MRII along with MR, if not, enter any other number:\n'))


fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize=(15, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()
xlab = r'$log_{10}(M_{dust}/M_{\odot})$'
ylab = r'$log_{10}(\Phi/(Mpc^{-3}dex^{-1}))$'   
ylim = [-6, -0.1]

if turn_on == 1:
    
    xlim = [4, 10.4]
    xticks = [4, 5, 6, 7, 8, 9, 10]
    savename = 'dust_mass_fn_MR_MRII.png'
    
else:

    xlim = [6, 10.4]
    xticks = [6, 7, 8, 9, 10]
    savename = 'dust_mass_fn_MR.png'


from obs_plots import DMF

for z in range(0, 9):

    file_req = files_(turn_on, z)
    
    if np.isscalar(file_req):
    
        xx, yy = dust_massfn(file_req, MR_vol)
        
        axs[z].plot(xx, np.log10(yy), color = 'black', lw = 2)
        
    else:
        xx, yy = dust_massfn(file_req[0], MR_vol)
        
        axs[z].plot(xx, np.log10(yy), color = 'black', lw = 2)
        
        xx, yy = dust_massfn(file_req[1], MRII_vol)
        
        axs[z].plot(xx, np.log10(yy), color = 'brown', lw = 2)
    
    
    DMF(axs[z], z)   #Plotting the observational data points
    
    if (z == 0) and (turn_on == 0):
        axs[z].text(6.5, -5, r'$z = {}$'.format(z), fontsize = 18)
    else:
        axs[z].text(9.5, -1, r'$z = {}$'.format(z), fontsize = 18)
    
    axs[z].set_xlim(xlim)
    axs[z].set_ylim(ylim)
    axs[z].set_xticks(xticks)
    axs[z].grid()
    axs[z].legend(frameon=False, fontsize = 16, markerscale=2, loc ='best', numpoints=1, handletextpad=0.005)
    for label in (axs[z].get_xticklabels() + axs[z].get_yticklabels()):
        label.set_fontsize(18)
    

fig.tight_layout()    
fig.subplots_adjust(bottom=0.09, left = 0.08, wspace=0, hspace=0)
fig.text(0.03, 0.5, ylab, va='center', rotation='vertical', fontsize=22)
fig.text(0.5, 0.04, xlab, va='center', fontsize=22)
plt.savefig(savename)
plt.show()
