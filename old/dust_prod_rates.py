import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
sys.path.append('/lustre/scratch/astro/ap629/my_lib/')
from essent import *

h = 0.673 #little h as used in the model


inp = int(input('Enter the appropriate number for the required dust production rate: \n0 - Cosmic dust rate versus redshift \n1 - Dust production rate versus Stellar Mass for different redshifts\n'))

if inp == 0:
    xlab = r'$z$'
    ylab = r'$log_{10}(\Phi_{DR}/(M_{\odot}yr^{-1}Mpc^{-3}))$'
    savename2 = 'dust_prod_rate_redshift'
elif inp == 1:
    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(\Phi_{DR}/(M_{\odot}yr^{-1}))$'
    savename2 = 'dust_prod_rate_Mstar'
else:
    print ('Not an applicable choice, retry....')
    sys.exit()  

norm_z = np.arange(0, 9, 1)
high_z = np.arange(9, 14, 1)
all_z = np.arange(0, 14, 1)

z_inp = int(input('Select redshift (z) range: \n0 - [0, 8]\n1 - [9, 13]\n2 - [0, 13]\n'))

if z_inp == 0:
    z_ = norm_z
    savename1 = 'z0_8'
    r = c = 3
    xlim = [6,11.5]
    ylim = [-6, 2.5]
    xticks = [6, 7, 8, 9, 10, 11, 12]
elif z_inp == 1:
    z_ = high_z
    savename1 = 'z9_13'
    r = c = 2
    xlim = [6,10.5]
    ylim = [-8, 1]
    xticks = [6, 7, 8, 9, 10]
elif z_inp == 2:
    z_ = all_z
    savename1 = 'z0_13'
    r = 5 
    c = 3
    xlim = [6,11.5]
    ylim = [-6, 2.5]
    xticks = [6, 7, 8, 9, 10, 11, 12]
else:
    print ('Not an applicable choice, retry....')
    sys.exit()  
 
turn_on = int(input('Input 1 to plot MRII, for MR any would do. \n'))

if turn_on == 1:
    savename3 = '_MRII.png'
else:
    savename3 = '_MR.png'
    
def files_(turn_on, z):

    if turn_on == 1:
        file_req = '/lustre/scratch/astro/ap629/Dust_output_17jan/MRII/hdf5/lgal_z{}_N*.hdf5'.format(z)
    else:
        file_req = '/lustre/scratch/astro/ap629/Dust_output_17jan/MR/hdf5/lgal_z{}_N*.hdf5'.format(z)
        #file_req = '/lustre/scratch/astro/ap629/sc558/data/MR/lgal_z{}_N*.hdf5'.format(z)
    return file_req
    
if inp == 0:

    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(9, 12), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    pAGB = pSNII = pSNIA = pGG = dest = pNET = pSFR = np.array([])
    for z in z_:
        
        file_req = files_(turn_on, z)

        
        Mstar = get_data1('StellarMass', file_req)*(1e10)/h
        Mdust = np.sum(get_data1('Dust_elements', file_req), axis = 1)
        SFR = get_data1('Sfr', file_req)
        sSFR = SFR/Mstar
        Type = get_data1('Type', file_req)
        
        vol = (480.279/h)**3
        
        dustrates = get_data1('DustRatesISM', file_req)
        ok = np.logical_and(sSFR > sSFR_cut(z), np.logical_and(Mdust > 0.0, Type == 0))
        #ok = np.logical_and(Mstar > 0, Mdust > 0)
        
        dustrates = dustrates[ok]
        SFR = SFR[ok]
        
        pAGB = np.append(pAGB, np.nansum(dustrates[:,0])/vol)
        pSNII = np.append(pSNII, np.nansum(dustrates[:,1])/vol)
        pSNIA = np.append(pSNIA, np.nansum(dustrates[:,2])/vol)
        pGG = np.append(pGG, np.nansum(dustrates[:,3])/vol)
        dest = np.append(dest, np.nansum(dustrates[:,4])/vol)
        pNET = np.append(pNET, np.nansum((np.nansum(dustrates[:, :4], axis = 1) - dustrates[:,4])/vol))
        pSFR = np.append(pSFR, np.nansum(SFR)/vol)
        
    axs.plot(z_, np.log10(pAGB), label = r'$AGB$', color='b', linewidth = 2)
    axs.plot(z_, np.log10(pSNII), label = r'$SNII$', color='r', linewidth = 2)
    axs.plot(z_, np.log10(pSNIA), label = r'$SNIA$', color='y', linewidth = 2)
    axs.plot(z_, np.log10(pGG), label = r'$GG$', color='g', linewidth = 2)
    axs.plot(z_, np.log10(dest), label = r'$DEST$', color='k', linewidth = 2)
    axs.plot(z_, np.log10(pNET), label = r'$NET$', color='orange', linewidth = 2)
    axs.plot(z_, np.log10(pSFR), label = r'$SFR$', color='cyan', linewidth = 2, ls = 'dashed')
    
    lgd = axs.legend(frameon=False, fontsize = 16, markerscale=2, loc = 3, numpoints=1, handletextpad=0.005)
    lgd.set_zorder(100)
    axs.grid()
    axs.set_xticks(z_)
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(18)

elif inp == 1:

    fig, axs = plt.subplots(nrows = r, ncols = c, figsize=(15, 13), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()
    
    for i, z in enumerate(z_):
        
        file_req = files_(turn_on, z)
        pAGB = pSNII = pSNIA = pGG = dest = np.array([])
        
        Mstar = get_data1('StellarMass', file_req)*(1e10)/h
        SFR = get_data1('Sfr', file_req)
        sSFR = SFR/Mstar
        Type = get_data1('Type', file_req)
        
        dustrates = get_data1('DustRatesISM', file_req)
        #ok = np.logical_and(Mstar > 0, Mstar > 0)
        ok = np.logical_and(sSFR > sSFR_cut(z), Type == 0)
        
        dustrates = dustrates[ok]
        SFR = SFR[ok]
        Mstar = Mstar[ok]
        
        xx, yy, yy_up, yy_low = get_median(Mstar, dustrates[:,0], n = 18)
        axs[i].plot(np.log10(xx), np.log10(yy), label = r'$AGB$', color='b', linewidth = 2)
        axs[i].plot(np.log10(xx), np.log10(yy_up), color='b', linewidth = 2, ls = 'dashed')
        axs[i].plot(np.log10(xx), np.log10(yy_low), color='b', linewidth = 2, ls = 'dashed')

        xx, yy, yy_up, yy_low = get_median(Mstar, dustrates[:,1], n = 18)
        axs[i].plot(np.log10(xx), np.log10(yy), label = r'$SNII$', color='r', linewidth = 2)
        axs[i].plot(np.log10(xx), np.log10(yy_up), color='r', linewidth = 2, ls = 'dashed')
        axs[i].plot(np.log10(xx), np.log10(yy_low), color='r', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_median(Mstar, dustrates[:,2], n = 18)
        axs[i].plot(np.log10(xx), np.log10(yy), label = r'$SNIA$', color='y', linewidth = 2)
        axs[i].plot(np.log10(xx), np.log10(yy_up), color='y', linewidth = 2, ls = 'dashed')
        axs[i].plot(np.log10(xx), np.log10(yy_low), color='y', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_median(Mstar, dustrates[:,3], n = 18)
        axs[i].plot(np.log10(xx), np.log10(yy), label = r'$Grain$ $Growth$', color='g', linewidth = 2)
        axs[i].plot(np.log10(xx), np.log10(yy_up), color='g', linewidth = 2, ls = 'dashed')
        axs[i].plot(np.log10(xx), np.log10(yy_low), color='g', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_median(Mstar, dustrates[:,4], n = 18)
        axs[i].plot(np.log10(xx), np.log10(yy), label = r'$Destruction$', color='k', linewidth = 2)
        axs[i].plot(np.log10(xx), np.log10(yy_up), color='k', linewidth = 2, ls = 'dashed')
        axs[i].plot(np.log10(xx), np.log10(yy_low), color='k', linewidth = 2, ls = 'dashed')
        
        axs[i].text(8.75, -5.5, r'$z = {}$'.format(z), fontsize = 18)
        
        
        if i == 7:
	        axs[i].set_xlim([9.5,11.95])
        axs[i].set_xticks(xticks)
        axs[i].set_xlim(xlim)
        axs[i].set_ylim(ylim)
        if z == z_[-1]:
            lgd = axs[i].legend(frameon=False, fontsize = 16, markerscale=2, loc = 4, numpoints=1, handletextpad=0.05)
            lgd.set_zorder(100)
        axs[i].grid()
        for label in (axs[i].get_xticklabels() + axs[i].get_yticklabels()):
            label.set_fontsize(18)

    if z_inp == 2:
        fig.delaxes(axs[-1])

fig.tight_layout()    
fig.subplots_adjust(bottom=0.1, left = 0.1, wspace=0, hspace=0)
fig.text(0.01, 0.5, ylab, va='center', rotation='vertical', fontsize=24)
fig.text(0.5, 0.04, xlab, va='center', fontsize=24)
plt.savefig(savename1+savename2+savename3)
plt.show()
