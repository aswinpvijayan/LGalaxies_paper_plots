import sys
import os
import warnings
#if not sys.warnoptions:
#    warnings.simplefilter("ignore")
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import get_
import seaborn as sns
sns.set_context("paper")

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
        file_req = '/lustre/scratch/astro/ap629/Dust_output_22May/MR/SA_output_*'
        #file_req = '/lustre/scratch/astro/ap629/sc558/data/MR/lgal_z{}_N*.hdf5'.format(z)
    return file_req

redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])

snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])


def dust_prod_rate_redshift(files, z):

    
    Mstar = (get_.get_var(files, 'StellarMass', snap)*1e10)/0.673
    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust1 = np.sum(Mdust1, axis = 1)  
    Mdust2 = np.sum(Mdust2, axis = 1) 
    Mdust = Mdust1 + Mdust2
    SFR = get_.get_var(files, 'Sfr', snap)
    sSFR = SFR/Mstar
    Type = get_.get_var(files, 'Type', snap)
    dustrates = get_.get_var(files, 'DustColdGasRates', snap)
    
    
    ok = np.logical_and(sSFR > get_.sSFR_cut(z), np.logical_and(Mdust > 0.0, Type == 0))
        
    dustrates = dustrates[ok]
    SFR = SFR[ok]
    Mstar = Mstar[ok]
    
    return dustrates, Mstar, SFR
    
if inp == 0:

    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(9, 12), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    pAGB = pSNII = pSNIA = pGG = dest = pNET = pSFR = np.array([])
    for z in z_:
        
        snap = snapnumMR[np.where(redshift == str(z))[0][0]]
        vol = (480.279/h)**3
        files = '/lustre/scratch/astro/ap629/Dust_output_8june/MR/SA_output_*'
        dustrates, Mstar, SFR = dust_prod_rate_redshift(files, z)
        
        dustrates = np.nansum(dustrates/vol, axis = 0)
        SFR = np.nansum(SFR/vol, axis = 0)        
        
        pAGB = np.append(pAGB, dustrates[0])
        pSNII = np.append(pSNII, dustrates[1])
        pSNIA = np.append(pSNIA, dustrates[2])
        pGG = np.append(pGG, dustrates[3])
        dest = np.append(dest, dustrates[4])
        pNET = np.append(pNET, np.nansum(dustrates[:4]) - dustrates[4])
        pSFR = np.append(pSFR, SFR)
        
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
        
        snap = snapnum[np.where(redshift == str(z))[0][0]]
        vol = (480.279/h)**3
        
        
        dustrates, Mstar, SFR = dust_prod_rate_redshift(files, z)
        
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
