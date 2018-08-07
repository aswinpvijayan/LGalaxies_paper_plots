#!/usr/bin/env python
import datetime
print ('Time at the start is: ', str(datetime.datetime.now()))

#This computes the dust mass density inside the simulation boxes.

import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
sys.path.append('/lustre/scratch/astro/ap629/my_lib/')
from essent import *

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
    
    xx, yy, yy_up, yy_low = get_median(thisx, thisy, n = 13)
    
    return xx, yy, yy_up, yy_low

h = 0.673 #little h as used in the model

#inp = int(input('Enter the appropriate number for the required phase space plot: \n0 - Stellar Mass versus Dust Mass \n1 - Metallicity versus Dust Mass\n2 - Dust fraction in cold gas versus Stellar Mass\n3 - Dust fraction in cold gas versus Metallicity\n'))

inp = int(sys.argv[1])

if inp == 0:
    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(M_{Dust}/(M_{\odot}))$'
    savename = 'Dust_Stellar_MR_MRII'
elif inp == 1:
    xlab = r'$12+log_{10}(O/H)$'
    ylab = r'$log_{10}(M_{Dust}/(M_{\odot}))$'
    savename = 'Dust_metal_MR_MRII'
elif inp == 2:
    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(M_{Dust}/(M_{cold\mathrm{ }gas}))$'
    savename = 'Dust_gas_ratio_stellar_mass_MR_MRII'
elif inp == 3:
    xlab = r'$12+log_{10}(O/H)$'
    ylab = r'$log_{10}(M_{Dust}/(M_{cold\mathrm{ }gas}))$'
    savename = 'Dust_gas_ratio_metal_MR_MRII'
else:
    print ('Not an applicable choice, retry....')
    sys.exit()  

"""
turn_on = int(input('Input 1 in order to plot MRII:\n'))

if turn_on ==  1:
    namex = '_MRII'
else:
    namex = ''
"""

Vars = np.array([['StellarMass', 'Dust_elements'], ['ColdGas_elements', 'Dust_elements'], ['ColdGas_elements', 'Dust_elements'], ['ColdGas_elements', 'Dust_elements']])

colours = ['blue', 'green', 'red', 'magenta', 'brown', 'orange', 'violet', 'cyan', 'indigo', 'black', 'yellow', 'blue', 'green', 'red']

fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize=(15, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

for z in range(0, 9):
    
    for turn_on in [0, 1]:
        
        if turn_on == 1:
            file_req = '/lustre/scratch/astro/ap629/Dust_output_17jan/MRII/hdf5/lgal_z{}_N*.hdf5'.format(z) 
            lab = 'MRII'
            num = 0
        else:
            file_req = '/research/astro/virgo/Aswin/data/MR/hdf5/lgal_z{}_N*.hdf5'.format(z)
            lab = 'MR'
            num = 1
        
        Mstar = get_data1('StellarMass', file_req)*(1e10)/h  #Now in units of Msun 
        SFR = get_data1('Sfr', file_req)    #Star formation rate, in units of Msun/yr
        sSFR = SFR/Mstar    #Specific star formation rate, in units of /yr
        del SFR
        Type = get_data1('Type', file_req)
        
        if inp == 0:
            
            from obs_plots import DM_obs
            
            Mdust = get_data1(Vars[inp][1], file_req)  #Array of mass of 11 elements in dust: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
            Mdust = np.sum(Mdust, axis = 1)  #Total dust mass
            
            xx, yy, yy_up, yy_low = outputs(Mstar, Mdust, sSFR, Type, z)
            
            xx = np.log10(xx)
            yy = np.log10(yy)
            yy_up = np.log10(yy_up)
            yy_low = np.log10(yy_low)
            
            xlim = [6.5,11.5]
            ylim = [2,9.95]
            xticks = [7, 8, 9, 10, 11, 12]
            axs[z].text(8.75, 9, r'$z = {}$'.format(z), fontsize = 18)
            
            #DM_obs(axs[z], z)   #Plotting the observational data points
            
            
        if inp == 1:
            
            from obs_plots import D_Met_obs
            
            Mcg = get_data1(Vars[inp][0], file_req)  #Array of mass of 11 elements present in the cold gas: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
            Mdust = get_data1(Vars[inp][1], file_req) #Array of mass of 11 elements in dust: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
            Mdust = np.sum(Mdust, axis = 1)  #Total dust mass
            met = Mcg[:,4]/(15.9994*Mcg[:,0])  #Metallicity = Number of Oxygen atoms/Number of Hydrogen atoms
            xx, yy, yy_up, yy_low = outputs(met, Mdust, sSFR, Type, z)
            
            xx = 12 + np.log10(xx)
            yy = np.log10(yy)
            yy_up = np.log10(yy_up)
            yy_low = np.log10(yy_low)
            
            xlim = [5.5,10]
            ylim = [2,9.95]
            xticks = [6, 7, 8, 9, 10]
            axs[z].text(5.75, 9, r'$z = {}$'.format(z), fontsize = 18)
            
            #D_Met_obs(axs[z], z)    #Plotting the observational data points
            
        if inp == 2:
            
            from obs_plots import DG_Mstar_obs
            
            Mcg = get_data1(Vars[inp][0], file_req) #Array of mass of 11 elements present in the cold gas: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
            Mcg = np.sum(Mcg, axis = 1) #Total cold gas mass
            Mdust = get_data1(Vars[inp][1], file_req) #Array of mass of 11 elements in dust: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
            Mdust = np.sum(Mdust, axis = 1) #Total dust mass
            Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
            xx, yy, yy_up, yy_low = outputs(Mstar, Mratio, sSFR, Type, z)
            
            xx = np.log10(xx)
            yy = np.log10(yy)
            yy_up = np.log10(yy_up)
            yy_low = np.log10(yy_low)
            
            xlim = [6.5,11.5]
            ylim = [-7,-0.1]
            xticks = [7, 8, 9, 10, 11]
            axs[z].text(6.75, -1, r'$z = {}$'.format(z), fontsize = 18)
            
            #DG_Mstar_obs(axs[z], z) #Plotting the observational data points
            
        if inp == 3:
            
            from obs_plots import DG_met_obs
            
            Mcg = get_data1(Vars[inp][0], file_req) #Array of mass of 11 elements present in the cold gas: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
            met = Mcg[:,4]/(15.9994*Mcg[:,0]) #Metallicity = Number of Oxygen atoms/Number of Hydrogen atoms
            Mcg = np.sum(Mcg, axis = 1) #Total cold gas mass
            Mdust = get_data1(Vars[inp][1], file_req) #Array of mass of 11 elements in dust: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
            Mdust = np.sum(Mdust, axis = 1)  #Total dust mass
            Mratio = Mdust/Mcg  #Dust-to-total cold gas mass ratio. Dimensionless.
            xx, yy, yy_up, yy_low = outputs(met, Mratio, sSFR, Type, z)
            
            xx = 12 + np.log10(xx)
            yy = np.log10(yy)
            yy_up = np.log10(yy_up)
            yy_low = np.log10(yy_low)
            
            xlim = [5.5,10]
            ylim = [-7,-0.1]
            xticks = [6, 7, 8, 9, 10]
            axs[z].text(5.75, -1, r'$z = {}$'.format(z), fontsize = 18)
            
            #DG_met_obs(axs[z], z)   #Plotting the observational data points
          
        
        axs[z].plot(xx, yy, lw = 2, color = colours[num], label = r'${}$'.format(lab))
        axs[z].plot(xx, yy_up, lw = 2, ls = 'dashed', color = colours[num])
        axs[z].plot(xx, yy_low, lw = 2, ls = 'dashed', color = colours[num])
        
        del xx, yy, yy_up, yy_low
        
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
plt.savefig(savename+'.png')
print ('End time is ', str(datetime.datetime.now()))
