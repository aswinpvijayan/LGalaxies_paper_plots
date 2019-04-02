#!/usr/bin/env python
import datetime
print ('Time at the start is: ', str(datetime.datetime.now()))

import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
sys.path.append('func_def/')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import get_
import seaborn as sns
sns.set_context("paper")
sns.set_style(style='white') 

        ###################################################################################
"""        
    inp = 0 to plot just MR, 1 for MRII and any higher number for plotting both MR and MRII
    (Figure 13 in paper)
"""        

inp = int(sys.argv[1])    

        ####################################################################################

redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])
sims = np.array(['MR', 'MRII'])

filesMR = '../Dust_output/MR/SA_output_*'
filesMRII = '../Dust_output/MRII/SA_output_*'
files = np.array([filesMR, filesMRII])

h = 0.673 #little h as used in the model

#Simulation volumes
MR_vol = (480.279/h)**3  #Millennium
MRII_vol = (96.0558/h)**3 #Millennium II
vol = np.array([MR_vol, MRII_vol])

def outputs(x, y, sSFR, Type, z):
    
    """
    #Function remove_ is for removing any nan's, inifnities and other non-physical values, can specify if the zero values in the array should not be removed, by giving the array number as argument. Returns the number density in a pixel and the median (in 13 bins in the x axis in logarithmic space) with the values of the 16th percentile and the 84th percentile.
    """
    
    out = get_.remove_(np.array([x, y, sSFR, Type]), np.array([3]))  
    out = out[(out[1] > 1e4) & (out[3] == 0)]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    this_sSFR = out[2]
    this_Type = out[3]
    
    del out, x, y, sSFR, Type
    
    return thisy


def dust_massfn(files, z, snapnum, i):

    """
    Outputs the dust mass function for input files of a particular redshift, for a given volume
    """

    add = sims[i] 
    snap = snapnum[i][np.where(redshift == str(z))[0][0]]
    Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
    Type = get_.get_var(files[i], 'Type', snap)
    Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    Mdust = Mdust1 + Mdust2
    
    SFR = get_.get_var(files[i], 'Sfr', snap)

    sSFR = SFR/Mstar
    
    ok = np.logical_and(Mdust > 0.0, Type == 0)
    #ok = np.logical_and(sSFR < get_.sSFR_cut(z), np.logical_and(Mdust > 0.0, Type == 0))
    out = outputs(Mstar[ok], Mdust[ok], sSFR[ok], Type[ok], z)
    dbins = np.arange(np.log10(min(out))-0.05, np.log10(max(out))+0.05, 0.2)
    bins, edges = np.histogram(np.log10(out), dbins)
    
    xx = (edges[1:]+edges[:-1])/2
    binsize = (xx[-1] - xx[0])/len(xx)
    yy = bins/(vol[i]*binsize)   
    
    return xx, yy 


fig, axs = plt.subplots(nrows = 2, ncols = 1, figsize=(6, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()
xlab = r'$\mathrm{log}_{10}(\mathrm{M}_{\mathrm{Dust}}/\mathrm{M}_{\odot})$'
ylab = r'$\mathrm{log}_{10}(\Phi/(\mathrm{M}_{\odot}\mathrm{yr}^{-1}\mathrm{Mpc}^{-3}))$'   
ylim = [-7, -0.1]

if inp == 0:
    
    xlim = [6, 10.4]
    xticks = [6, 7, 8, 9, 10]
    savename = 'dust_mass_fn_MR.pdf'
    print ('Plotting MR dust mass function')

elif inp == 1:
    
    xlim = [4, 10.4]
    xticks = [4, 5, 6, 7, 8, 9, 10]
    savename = 'dust_mass_fn_MRII.pdf'
    print ('Plotting MRII dust mass function')

else:
    
    xlim = [5, 10.5]
    xticks = [5, 6, 7, 8, 9, 10]
    savename = 'dust_mass_fn_MR_MRII.pdf'
    print ('Plotting MR (black) and MRII (brown) dust mass function')    

from obs_plots import DMF

for z in range(0, 2):

    if inp in [0, 1]:
    
        xx, yy = dust_massfn(files, z, snaps, inp)
        
        if z == 8:
            axs[z].plot(xx, np.log10(yy), color = 'black', lw = 2, legend = r'${}$'.format(sims[inp]))
        else:
            axs[z].plot(xx, np.log10(yy), color = 'black', lw = 2)
        
    else:

        xx, yy = dust_massfn(files, z, snaps, 0)
        
        if z == 8:
            axs[z].plot(xx, np.log10(yy), color = 'black', lw = 2, label = r'$\mathrm{MR}$')
        else:
            axs[z].plot(xx, np.log10(yy), color = 'black', lw = 2)
        
        xx, yy = dust_massfn(files, z, snaps, 1)
        
        if z == 8:
            axs[z].plot(xx, np.log10(yy), color = 'brown', lw = 2, label = r'$\mathrm{MRII}$')
        else:
            axs[z].plot(xx, np.log10(yy), color = 'brown', lw = 2)
    
    DMF(axs[z], z)   #Plotting the observational data points
    
    axs[z].set_xlim(xlim)
    axs[z].set_ylim(ylim)
    axs[z].set_xticks(xticks)
    axs[z].grid()
    axs[z].legend(frameon=False, fontsize = 16, markerscale=2, loc ='best', numpoints=1, handletextpad=0.005)
    for label in (axs[z].get_xticklabels() + axs[z].get_yticklabels()):
        label.set_fontsize(18)
    

fig.tight_layout()    
fig.subplots_adjust(bottom=0.09, left = 0.15, wspace=0, hspace=0)
fig.text(0.02, 0.51, ylab, va='center', rotation='vertical', fontsize=22)
fig.text(0.35, 0.04, xlab, va='center', fontsize=22)
plt.savefig(savename)

print ('End time is ', str(datetime.datetime.now()))
