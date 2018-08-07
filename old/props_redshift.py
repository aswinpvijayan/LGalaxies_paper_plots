#!/usr/bin/env python

#This computes the dust mass density inside the simulation boxes.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
sys.path.append('/lustre/scratch/astro/ap629/my_lib/')
from essent import *

h = 0.673 #little h as used in the model

MR_vol = (480.279/h)**3
MRII_vol = (96.0558/h)**3

def files_(z):

    file_req = ['/lustre/scratch/astro/ap629/Dust_output_17jan/MR/hdf5/lgal_z{}_N*.hdf5'.format(z), '/lustre/scratch/astro/ap629/Dust_output_17jan/MRII/hdf5/lgal_z{}_N*.hdf5'.format(z)] 

    return file_req

def dust_massdensity(z_, sim):

    """
    Outputs the dust production density for input files of a particular redshift, for a given box
    """
    
    den = np.ones(len(z_))*np.nan
    for i, z in enumerate(z_):
        if sim == 'MR':
            vol = MR_vol
            file_req = files_(z)[0]
        elif sim == 'MRII':
            vol = MR_vol
            file_req = files_(z)[1]
        Mstar = get_data1('StellarMass', file_req)*(1e10)/h  #Now in units of Msun 
        SFR = get_data1('Sfr', file_req)    #Star formation rate, in units of Msun/yr
        sSFR = SFR/Mstar    #Specific star formation rate, in units of /yr
        Type = get_data1('Type', file_req)
        Mdust = get_data1('Dust_elements', file_req)
        Mdust = np.sum(Mdust, axis = 1)
        out = remove_(np.array([Mdust, sSFR, Type]), np.array([2]))
        Mdust = out[0] 
        sSFR = out[1]
        Type = out[2]
        ok = np.logical_and(sSFR > sSFR_cut(z), Type == 0)
        den[i] = np.log10(np.sum(Mdust[ok]/vol))  
    
    return den


xlab = r'$z$'
ylab = r'$log_{10}(M_{dust}/(M_{\odot}Mpc^{-3}))$'

turn_on = int(input('Input 1 in order to include MRII:\n'))
if turn_on == 1:
    name_x = '_MRII'
else:
    name_x = ''

norm_z = np.arange(0, 9, 1)
high_z = np.arange(9, 14, 1)
all_z = np.arange(0, 14, 1)

z_inp = int(input('Select redshift (z) range: \n0 - [0, 8]\n1 - [9, 13]\n2 - [0, 13]\n'))

if z_inp == 0:
    z_ = norm_z
    savename = 'z0_8'
elif z_inp == 1:
    z_ = high_z
    savename = 'z9_13'
elif z_inp == 2:
    z_ = all_z
    savename = 'z0_13'
else:
    print ('Not an applicable choice, retry....')
    sys.exit()  


fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(11, 12), sharex=True, sharey=True, facecolor='w', edgecolor='k')


if turn_on == 1:
    
    axs.plot(z_, dust_massdensity(z_, 'MR'), color = 'black', lw = 2)
    axs.plot(z_, dust_massdensity(z_, 'MRII'), color = 'brown', lw = 2)

else:

    axs.plot(z_, dust_massdensity(z_, 'MR'), color = 'black', lw = 2)    

axs.set_xticks(z_)
axs.grid()
for label in (axs.get_xticklabels() + axs.get_yticklabels()):
    label.set_fontsize(18)

fig.tight_layout()    
fig.subplots_adjust(bottom=0.09, left = 0.1, wspace=0, hspace=0)
fig.text(0.01, 0.5, ylab, va='center', rotation='vertical', fontsize=22)
fig.text(0.5, 0.04, xlab, va='center', fontsize=26)
plt.savefig(savename+'dust_prod_density_MR'+name_x+'.png')
plt.show()

