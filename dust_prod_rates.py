import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
sys.path.append('func_def/')
import matplotlib as mpl
mpl.use('Agg') 
mpl.rcParams['text.usetex'] = False  
import matplotlib.pyplot as plt
import get_
import seaborn as sns
sns.set_context("paper")
sns.set_style(style='white') 

            ##########################################################################
"""
    User inputs for producing the preferred plots: 
    0 - Cosmic dust rate versus redshift    (Figure 6 in paper)
    1 - Dust production rate versus Stellar Mass for different redshifts    (Figure 8 in paper)
    1 - Dust production rate versus Metallicity for different redshifts
    
    inp = user input for preferred plot
    
    z_inp = For plotting different redshift range. 
    0 - [0, 8]  (Figure 6, 8 in paper)
    1 - [9, 13] (Figure 13c in paper)
    2 - [0, 13]
            
    i = 0 to plot just MR, 1 for MRII and any higher number for plotting both MR and MRII        
"""            

h = 0.673 #little h as used in the model

filesMR = '../Dust_output/MR/SA_output_*'
filesMRII = '../Dust_output/MRII/SA_output_*'
files = np.array([filesMR, filesMRII])

redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])

snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])
snaps = np.array([snapnumMR, snapnumMRII])
sims = np.array(['MR', 'MRII'])

vol = np.array([(480.279/h)**3, (96./h)**3])

inp = int(sys.argv[1])
z_inp = int(sys.argv[2])
i = int(sys.argv[3])

            ##########################################################################

norm_z = np.arange(0, 9, 1)
high_z = np.arange(9, 13, 1)
all_z = np.arange(0, 13, 1)

if i == 1:
    savename3 = '_MRII'
elif i == 0:
    savename3 = '_MR'
else:
    savename3 = 'MR_MRII'

if inp == 0:
    xlab = r'$\mathrm{z}$'
    ylab = r'$\mathrm{log}_{10}(\Phi_{\mathrm{DR}}/(M_{\odot}\mathrm{yr}^{-1}\mathrm{Mpc}^{-3}))$'
    savename2 = 'dust_prod_rate_redshift'
elif inp == 1:
    xlab = r'$\mathrm{log}_{10}(M_{*}/M_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(\Phi_{\mathrm{DR}}/(M_{\odot}\mathrm{yr}^{-1}))$'
    savename2 = 'dust_prod_rate_Mstar'
elif inp == 2:
    xlab = r'$12+\mathrm{log}_{10}(\mathrm{O/H})$'
    ylab = r'$\mathrm{log}_{10}(\Phi_{\mathrm{DR}}/(M_{\odot}\mathrm{yr}^{-1}))$'
    savename2 = 'dust_prod_rate_Z'
else:
    print ('Not an applicable choice, retry....')
    sys.exit()  

if z_inp == 0:
    z_ = norm_z
    savename1 = 'z0_8'
    r = c = 3
    if inp == 0 or inp == 1:
        xlim = [7.5,11.5]
        xticks = [8, 9, 10, 11]
        ylim = [-6, 1.5]
    else:
        xlim = [7.3,10.7]
        xticks = [8, 9, 10]
        ylim = [-8, 1.5]
    
elif z_inp == 1:
    z_ = high_z
    savename1 = 'z9_12'
    r = c = 2
    if inp == 0 or inp == 1:
        xlim = [7.5,11.8]
        xticks = [8, 9, 10, 11]
        ylim = [-6, 2.5]
    else:
        xlim = [7.3,10.7]
        xticks = [8, 9, 10]
        ylim = [-8, 1.5]

elif z_inp == 2:
    z_ = all_z
    savename1 = 'z0_12'
    r = 5 
    c = 3
    if inp == 0 or inp == 1:
        xlim = [6,11.8]
        xticks = [6, 7, 8, 9, 10, 11]
        ylim = [-6, 1.5]
    else:
        xlim = [7.3,10.7]
        xticks = [8, 9, 10]
        ylim = [-8, 1.5]
else:
    print ('Not an applicable choice, retry....')
    sys.exit()  
 
def dust_prod_rate_redshift(files, z, snap, i):

    
    Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
    Type = get_.get_var(files[i], 'Type', snap)
    Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
    Mdust = Mdust1 + Mdust2
    SFR = get_.get_var(files[i], 'Sfr', snap)
    sSFR = SFR/Mstar
    dustrates = get_.get_var(files[i], 'DustColdGasRates', snap)
    Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)
    Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)
    Mcg = Mcg1 + Mcg2
    Z = (Mcg[:,4]-Mdust[:,4])/(15.9994*(Mcg[:,0]))  #O/H ratio
    #Z = np.nansum(Mcg[:,2:], axis = 1)/np.nansum(Mcg, axis = 1)
    Mdust = np.nansum(Mdust, axis = 1)
    Mcg = np.nansum(Mcg, axis = 1)
    ok = np.logical_and(sSFR > get_.sSFR_cut(z), np.logical_and(Mdust > 0.0, Type == 0))
        
    dustrates = dustrates[ok]
    SFR = SFR[ok]
    Mstar = Mstar[ok]
    Z = Z[ok]
    
    return dustrates, Mstar, SFR, Z
    
if inp == 0:
    if z_inp == 1:
        fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(6, 6), sharex=True, sharey=True, facecolor='w', edgecolor='k')
        bottom = 0.125
        left = 0.17
    else:
        fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(9, 12), sharex=True, sharey=True, facecolor='w', edgecolor='k')
        bottom = 0.1
        left = 0.11
    
    pAGB = pSNII = pSNIA = pGG = dest = pNET = pSFR = np.array([])
    for z in z_:
        snap = snaps[i][np.where(redshift == str(z))[0][0]]
        dustrates, Mstar, SFR, Z = dust_prod_rate_redshift(files, z, snap, i)
        ok = np.where(Mstar >= 10**8.9)
        dustrates = np.nansum(dustrates[ok]/vol[i], axis = 0)
        SFR = np.nansum(SFR[ok]/vol[i], axis = 0)        
        
        pAGB = np.append(pAGB, dustrates[0])
        pSNII = np.append(pSNII, dustrates[1])
        pSNIA = np.append(pSNIA, dustrates[2])
        pGG = np.append(pGG, dustrates[3])
        dest = np.append(dest, dustrates[4])
        pNET = np.append(pNET, np.nansum(dustrates[:4]) - dustrates[4])
        pSFR = np.append(pSFR, SFR)
        
    axs.plot(z_, np.log10(pAGB), label = r'$\mathrm{AGB}$', color='b', linewidth = 2)
    axs.plot(z_, np.log10(pSNII), label = r'$\mathrm{SNII}$', color='r', linewidth = 2)
    axs.plot(z_, np.log10(pSNIA), label = r'$\mathrm{SNIA}$', color='y', linewidth = 2)
    axs.plot(z_, np.log10(pGG), label = r'$\mathrm{GG}$', color='g', linewidth = 2)
    axs.plot(z_, np.log10(dest), label = r'$\mathrm{DEST}$', color='k', linewidth = 2)
    axs.plot(z_, np.log10(pNET), label = r'$\mathrm{NET}$', color='orange', linewidth = 3)
    axs.plot(z_, np.log10(pSFR), label = r'$\mathrm{SFR}$', color='cyan', linewidth = 2, ls = 'dashed')
    
    lgd = axs.legend(frameon=False, fontsize = 16, markerscale=2, loc = 3, numpoints=1, handletextpad=0.005)
    lgd.set_zorder(100)
    axs.grid()
    axs.set_xticks(z_)
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(18)

elif inp == 1:

    fig, axs = plt.subplots(nrows = 1, ncols = 3, figsize=(15, 6), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()
    bottom=0.13 
    left = 0.08
    for j, z in enumerate([0, 7, 8]):
        
        if i in [0,1]:
            snap = snaps[i][np.where(redshift == str(z))[0][0]]
            dustrates, Mstar, SFR, Z = dust_prod_rate_redshift(files, z, snap, i)
        else:
            snap = snaps[0][np.where(redshift == str(z))[0][0]]
            dustrates, Mstar, SFR, Z = dust_prod_rate_redshift(files, z, snap, 0)
            ok = np.where(Mstar > 10**9.0)
            dustrates = dustrates[ok]
            Mstar = Mstar[ok]
            SFR = SFR[ok]
            Z = Z[ok]
            snap = snaps[1][np.where(redshift == str(z))[0][0]]
            tmp1, tmp2, tmp3, tmp4 = dust_prod_rate_redshift(files, z, snap, 1)
            ok = np.logical_and(tmp2 > 10**7.0, tmp2 < 10**9.0)
            
            dustrates = np.append(dustrates, tmp1[ok], axis = 0)
            Mstar = np.append(Mstar, tmp2[ok], axis = 0)
            SFR = np.append(SFR, tmp3[ok], axis = 0)
            Z = np.append(Z, tmp4[ok], axis = 0)
            
        
        x = Mstar
        
        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,0], n = 12)
        axs[j].plot(np.log10(xx), np.log10(yy), label = r'$\mathrm{AGB}$', color='b', linewidth = 2)
        axs[j].plot(np.log10(xx), np.log10(yy_up), color='b', linewidth = 2, ls = 'dashed')
        axs[j].plot(np.log10(xx), np.log10(yy_low), color='b', linewidth = 2, ls = 'dashed')

        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,1], n = 12)
        axs[j].plot(np.log10(xx), np.log10(yy), label = r'$\mathrm{SNII}$', color='r', linewidth = 2)
        axs[j].plot(np.log10(xx), np.log10(yy_up), color='r', linewidth = 2, ls = 'dashed')
        axs[j].plot(np.log10(xx), np.log10(yy_low), color='r', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,2], n = 12)
        axs[j].plot(np.log10(xx), np.log10(yy), label = r'$\mathrm{SNIA}$', color='y', linewidth = 2)
        axs[j].plot(np.log10(xx), np.log10(yy_up), color='y', linewidth = 2, ls = 'dashed')
        axs[j].plot(np.log10(xx), np.log10(yy_low), color='y', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,3], n = 12)
        axs[j].plot(np.log10(xx), np.log10(yy), label = r'$\mathrm{Grain}$ $\mathrm{Growth}$', color='g', linewidth = 2)
        axs[j].plot(np.log10(xx), np.log10(yy_up), color='g', linewidth = 2, ls = 'dashed')
        axs[j].plot(np.log10(xx), np.log10(yy_low), color='g', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,4], n = 12)
        axs[j].plot(np.log10(xx), np.log10(yy), label = r'$\mathrm{Destruction}$', color='k', linewidth = 2)
        axs[j].plot(np.log10(xx), np.log10(yy_up), color='k', linewidth = 2, ls = 'dashed')
        axs[j].plot(np.log10(xx), np.log10(yy_low), color='k', linewidth = 2, ls = 'dashed')
        
        #xx, yy, yy_up, yy_low = get_.get_median(x, np.nansum(dustrates[:,0:4], axis = 1) - dustrates[:,4], n = 13)
        #axs[j].plot(np.log10(xx), np.log10(yy), label = r'$Net$', color='brown', linewidth = 2)
        #axs[j].plot(np.log10(xx), np.log10(yy_up), color='brown', linewidth = 2, ls = 'dashed')
        #axs[j].plot(np.log10(xx), np.log10(yy_low), color='brown', linewidth = 2, ls = 'dashed')
        
        axs[j].text(9.5, 0.5, r'$z = {}$'.format(z), fontsize = 18)
        
        
        axs[j].set_xticks(xticks)
        axs[j].set_xlim(xlim)
        axs[j].set_ylim(ylim)
        if z == z_[-1]:
            lgd = axs[j].legend(frameon=False, fontsize = 16, markerscale=2, loc = 0, numpoints=1, handletextpad=0.05)
            lgd.set_zorder(100)
        axs[j].grid()
        for label in (axs[j].get_xticklabels() + axs[j].get_yticklabels()):
            label.set_fontsize(18)

    if z_inp == 2:
        fig.delaxes(axs[-1])


elif inp == 2:

    fig, axs = plt.subplots(nrows = r, ncols = c, figsize=(15, 13), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()
    bottom=0.13 
    left = 0.08
    for j, z in enumerate(z_):
        
        if i in [0,1]:
            snap = snaps[i][np.where(redshift == str(z))[0][0]]
            dustrates, Mstar, SFR, Z = dust_prod_rate_redshift(files, z, snap, i)
        else:
            snap = snaps[0][np.where(redshift == str(z))[0][0]]
            dustrates, Mstar, SFR, Z = dust_prod_rate_redshift(files, z, snap, 0)
            ok0 = np.logical_and(Mstar >= 10**8.9, 12. + np.log10(Z) > 6.)
            snap = snaps[1][np.where(redshift == str(z))[0][0]]
            tmp1, tmp2, tmp3, tmp4 = dust_prod_rate_redshift(files, z, snap, 1)
            ok1 = np.logical_and(tmp2 <= 10**8.85, np.logical_and(tmp2 > 10**5.95, 12. + np.log10(tmp4) > 6.))
            
            dustrates = np.append(dustrates[ok0], tmp1[ok1], axis = 0)
            Mstar = np.append(Mstar[ok0], tmp2[ok1], axis = 0)
            SFR = np.append(SFR[ok0], tmp3[ok1], axis = 0)
            Z = np.append(Z[ok0], tmp4[ok1], axis = 0)
            
        
        x = Z
        
        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,0], n = 13)
        axs[j].plot(12. + np.log10(xx), np.log10(yy), label = r'$\mathrm{AGB}$', color='b', linewidth = 2)
        axs[j].plot(12. + np.log10(xx), np.log10(yy_up), color='b', linewidth = 2, ls = 'dashed')
        axs[j].plot(12. + np.log10(xx), np.log10(yy_low), color='b', linewidth = 2, ls = 'dashed')

        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,1], n = 13)
        axs[j].plot(12. + np.log10(xx), np.log10(yy), label = r'$\mathrm{SNII}$', color='r', linewidth = 2)
        axs[j].plot(12. + np.log10(xx), np.log10(yy_up), color='r', linewidth = 2, ls = 'dashed')
        axs[j].plot(12. + np.log10(xx), np.log10(yy_low), color='r', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,2], n = 13)
        axs[j].plot(12. + np.log10(xx), np.log10(yy), label = r'$\mathrm{SNIA}$', color='y', linewidth = 2)
        axs[j].plot(12. + np.log10(xx), np.log10(yy_up), color='y', linewidth = 2, ls = 'dashed')
        axs[j].plot(12. + np.log10(xx), np.log10(yy_low), color='y', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,3], n = 13)
        axs[j].plot(12. + np.log10(xx), np.log10(yy), label = r'$\mathrm{Grain}$ $\mathrm{Growth}$', color='g', linewidth = 2)
        axs[j].plot(12. + np.log10(xx), np.log10(yy_up), color='g', linewidth = 2, ls = 'dashed')
        axs[j].plot(12. + np.log10(xx), np.log10(yy_low), color='g', linewidth = 2, ls = 'dashed')
        
        xx, yy, yy_up, yy_low = get_.get_median(x, dustrates[:,4], n = 13)
        axs[j].plot(12. + np.log10(xx), np.log10(yy), label = r'$\mathrm{Destruction}$', color='k', linewidth = 2)
        axs[j].plot(12. + np.log10(xx), np.log10(yy_up), color='k', linewidth = 2, ls = 'dashed')
        axs[j].plot(12. + np.log10(xx), np.log10(yy_low), color='k', linewidth = 2, ls = 'dashed')
        
        #xx, yy, yy_up, yy_low = get_.get_median(x, np.nansum(dustrates[:,0:4], axis = 1) - dustrates[:,4], n = 13)
        #axs[j].plot(12. + np.log10(xx), np.log10(yy), label = r'$Net$', color='brown', linewidth = 2)
        #axs[j].plot(12. + np.log10(xx), np.log10(yy_up), color='brown', linewidth = 2, ls = 'dashed')
        #axs[j].plot(12. + np.log10(xx), np.log10(yy_low), color='brown', linewidth = 2, ls = 'dashed')
        
        axs[j].text(7.5, 1.2, r'$z = {}$'.format(z), fontsize = 18)
        
        
        axs[j].set_xticks(xticks)
        axs[j].set_xlim(xlim)
        axs[j].set_ylim(ylim)
        if z == z_[-1]:
            lgd = axs[j].legend(frameon=False, fontsize = 16, markerscale=2, loc = 4, numpoints=1, handletextpad=0.05)
            lgd.set_zorder(100)
        axs[j].grid()
        for label in (axs[j].get_xticklabels() + axs[j].get_yticklabels()):
            label.set_fontsize(18)

    if z_inp == 2:
        fig.delaxes(axs[-1])


fig.tight_layout()    
fig.subplots_adjust(bottom = bottom, left = left, wspace=0, hspace=0)
fig.text(0.01, 0.52, ylab, va='center', rotation='vertical', fontsize=24)
fig.text(0.49, 0.05, xlab, va='center', fontsize=26)
plt.savefig(savename1+savename2+savename3+'_full.pdf')
plt.close()
