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
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
import get_
import mbb
import gc
import seaborn as sns
sns.set_context("paper")
from astropy.cosmology import Planck13
from mpl_toolkits.axes_grid.inset_locator import inset_axes

        ###################################################################################
"""
    User inputs for producing the preferred plots:
    0. Dust mass vs stellar mass plot for z = 0-9   
    1. Dust mass vs Metallciity plot for z = 0-9
    2. Dust-to-gas ratio (DGR) vs stellar mass for z = 0-9
    3. DGR vs Metallicity for z = 0-9
    4. Dust-to-metal (DTM) ratio vs Stellar mass for z = 0-9  (Figure 3 in paper)
    5. DTM ratio vs Metallicity for z = 0-9
    6. Accretion timescale vs stellar mass plot for z = 0-9
    7. Accretion timescale vs metallicity plot for z = 0-9

    inp = user input for preferred plot

    i = 0 to plot just MR, 1 for MRII and any higher number for plotting both MR and MRII
"""
filesMR = '../Rob_dust_output/MR/SA_output_*'
filesMRII = '../Rob_dust_output/MRII/SA_output_*'
files = np.array([filesMR, filesMRII])

inp, i = int(sys.argv[1]), int(sys.argv[2])

        ###################################################################################

h = 0.673
redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])
sims = np.array(['MR', 'MRII'])

MR_vol = (480.279/h)**3  #Millennium
MRII_vol = (96.0558/h)**3 #Millennium II

fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize=(15, 13), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

if inp == 0:

    xlab = r'$\mathrm{log}_{10}(M_{*}/M_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(M_{\mathrm{Dust}}/M_{\odot})$'
    savename = 'Dust_Stellar_'

    from obs_plots import DM_obs
    from Mstar_Mdust import plot_Mstar_Mdust_age

    for z in range(0, 9):
    #    for k in range(0, len(titles)):
        DM_obs(axs[z], z)   #Plotting the observational data points
        if i <= 1:
            add, p, den = plot_Mstar_Mdust_age(files, z, axs[z], snaps, i, False)
        else:
            print ('Plotting dust mass vs stellar mass for z = {}'.format(z))
            add1, p, den = plot_Mstar_Mdust_age(files, z, axs[z], snaps, 0, True)
            add2, p, den = plot_Mstar_Mdust_age(files, z, axs[z], snaps, 1, True)
            add = add1+'_'+add2

        cbaxes = inset_axes(axs[z], width="95%", height="2%", loc=9)
        fig.colorbar(p, cax=cbaxes, orientation='horizontal')

        cbaxes.set_xticklabels([max(0.01, np.round(x, 2)) for x in cbaxes.get_xticks()*max(den)], fontsize=13)
        cbaxes.set_zorder(1)
        axs[z].text(8.45, 2.5, r'$z = {}$'.format(z), fontsize = 19)

elif inp == 1:

    xlab = r'$12+\mathrm{log}_{10}(\mathrm{O/H})$'
    ylab = r'$\mathrm{log}_{10}(M_{\mathrm{Dust}}/M_{\odot})$'
    savename = 'Dust_metal_'

    from obs_plots import D_Met_obs
    from Met_Mdust import plot_O_H_vs_Dust_age

    for z in range(0, 9):
        D_Met_obs(axs[z], z)    #Plotting the observational data points
        if i <= 1:
            add, p, den = plot_O_H_vs_Dust_age(files, z, axs[z], snaps, i, False)
        else:
            add, p, den = plot_O_H_vs_Dust_age(files, z, axs[z], snaps, i, True)

        cbaxes = inset_axes(axs[z], width="95%", height="2%", loc=9)
        fig.colorbar(p, cax=cbaxes, orientation='horizontal')

        cbaxes.set_xticklabels([max(0.01, np.round(x, 2)) for x in cbaxes.get_xticks()*max(den)], fontsize=13)
        cbaxes.set_zorder(1)
        axs[z].text(7.75, 9, r'$z = {}$'.format(z), fontsize = 18)

elif inp == 2:

    xlab = r'$\mathrm{log}_{10}(M_{*}/M_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(M_{\mathrm{Dust}}/M_{\mathrm{Cold gas}})$'
    savename = 'DGR_stell_'

    from obs_plots import DG_Mstar_obs
    from DGR_Mstar import plot_DGR_Mstar_age

    for z in range(0, 9):
        DG_Mstar_obs(axs[z], z) #Plotting the observational data points
        if i <= 1:
            add, p, den = plot_DGR_Mstar_age(files, z, axs[z], snaps, i, False)
        else:
            add1, p, den = plot_DGR_Mstar_age(files, z, axs[z], snaps, 0, True)
            add2, p, den = plot_DGR_Mstar_age(files, z, axs[z], snaps, 1, True)
            add = add1+'_'+add2

        cbaxes = inset_axes(axs[z], width="95%", height="2%", loc=9)
        fig.colorbar(p, cax=cbaxes, orientation='horizontal')

        cbaxes.set_xticklabels([max(0.01, np.round(x, 2)) for x in cbaxes.get_xticks()*max(den)], fontsize=13)
        cbaxes.set_zorder(1)
        axs[z].text(8.75, -1, r'$z = {}$'.format(z), fontsize = 18)

elif inp == 3:

    xlab = r'$12+\mathrm{log}_{10}(O/H)$'
    ylab = r'$\mathrm{log}_{10}(M_{\mathrm{Dust}}/M_{\mathrm{Cold gas}})$'
    savename = 'DGR_met_ratio_'

    from obs_plots import DG_met_obs
    from DGR_Met import plot_O_H_vs_DGR_age

    for z in range(0, 9):
        DG_met_obs(axs[z], z)   #Plotting the observational data points
        if i <= 1:
            add, p, den = plot_O_H_vs_DGR_age(files, z, axs[z], snaps, i, False)
        else:
            add, p, den = plot_O_H_vs_DGR_age(files, z, axs[z], snaps, i, True)

        cbaxes = inset_axes(axs[z], width="95%", height="2%", loc=9)
        fig.colorbar(p, cax=cbaxes, orientation='horizontal')

        cbaxes.set_xticklabels([max(0.01, np.round(x, 2)) for x in cbaxes.get_xticks()*max(den)], fontsize=13)
        cbaxes.set_zorder(1)
        axs[z].text(10, -5, r'$z = {}$'.format(z), fontsize = 18)

elif inp == 4:

    xlab = r'$\mathrm{log}_{10}(M_{*}/M_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(M_{\mathrm{Dust}}/M_{\mathrm{Metal}})$'
    savename = 'DTM_stellar_'

    from obs_plots import DTM_stell
    from DTM_Mstar import plot_Mstar_DTM_age

    for z in range(0, 9):
        DTM_stell(axs[z], z)
        if i <= 1:
            add, p, den = plot_Mstar_DTM_age(files, z, axs[z], snaps, i, False)
        else:
            add1, p, den = plot_Mstar_DTM_age(files, z, axs[z], snaps, 0, True)
            add2, p, den = plot_Mstar_DTM_age(files, z, axs[z], snaps, 1, True)
            add = add1+'_'+add2

        cbaxes = inset_axes(axs[z], width="93%", height="3%", loc=9)
        fig.colorbar(p, cax=cbaxes, orientation='horizontal')
        cbaxes.set_xticklabels([max(0.01, np.round(x, 2)) for x in cbaxes.get_xticks()*max(den)], fontsize=13)
        cbaxes.set_zorder(1)
        axs[z].text(11, -1.5, r'$z = {}$'.format(z), fontsize = 18)

elif inp == 5:

    xlab = r'$12+\mathrm{log}_{10}(O/H)$'
    ylab = r'$\mathrm{log}_{10}(M_{\mathrm{Dust}}/M_{\mathrm{Metal}})$'
    savename = 'DTM_met_'

    from obs_plots import DTM_oxy
    from DTM_Met import plot_O_H_vs_DTM_age

    for z in range(0, 1):
        DTM_oxy(axs[z], z)
        if i <= 1:
            add, p, den = plot_O_H_vs_DTM_age(files, z, axs[z], snaps, i, False)
        else:
            add, p, den = plot_O_H_vs_DTM_age(files, z, axs[z], snaps, i, True)

        cbaxes = inset_axes(axs[z], width="95%", height="2%", loc=9)
        fig.colorbar(p, cax=cbaxes, orientation='horizontal')

        cbaxes.set_xticklabels([max(0.01, np.round(x, 2)) for x in cbaxes.get_xticks()*max(den)], fontsize=13)
        cbaxes.set_zorder(1)
        axs[z].text(10, -1, r'$z = {}$'.format(z), fontsize = 18)


elif inp == 6:

    xlab = r'$\mathrm{log}_{10}(M_{*}/M_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(\tau_{\mathrm{acc}})$'
    savename = 'Mstar_tacc_'

    from tacc_Mstar import plot_tacc_Mstar_age

    for j,z in enumerate(range(0, 9)):

        uni_age = Planck13.age(z).value*1e9
        axs[j].axhline(y = np.log10(uni_age), ls = '-.',label = r'$\mathrm{Universe}$ $\mathrm{age(z)}$')
        if i <= 1:
            add, p, den = plot_tacc_Mstar(files, z, axs[j], snaps, i, False)
        else:
            add1, p, den = plot_tacc_Mstar(files, z, axs[j], snaps, 0, True)
            add2, p, den = plot_tacc_Mstar(files, z, axs[j], snaps, 1, True)
            add = add1+'_'+add2

        cbaxes = inset_axes(axs[z], width="95%", height="2%", loc=9)
        fig.colorbar(p, cax=cbaxes, orientation='horizontal')

        cbaxes.set_xticklabels([max(0.01, np.round(x, 2)) for x in cbaxes.get_xticks()*max(den)], fontsize=13)
        cbaxes.set_zorder(1)
        axs[z].text(8.2, 5.7, r'$z = {}$'.format(z), fontsize = 18)
        if z != 0:
            axs[z].legend().set_visible(False)

elif inp == 7:

    xlab = r'$12+\mathrm{log}_{10}(\mathrm{O/H})$'
    ylab = r'$\mathrm{log}_{10}(\tau_{\mathrm{acc}})$'
    savename = 'Met_tacc_'

    from tacc_Met import plot_tacc_Met_age

    for j,z in enumerate(range(0, 9)):

        uni_age = Planck13.age(z).value*1e9
        axs[j].axhline(y = np.log10(uni_age), ls = '-.',label = r'$\mathrm{Universe}$ $\mathrm{age(z)}$')
        if i <= 1:
            add, p, den = plot_tacc_Met_age(files, z, axs[z], snaps, i, False)
        else:
            add, p, den = plot_tacc_Met_age(files, z, axs[z], snaps, i, True)

        cbaxes = inset_axes(axs[z], width="95%", height="2%", loc=9)
        fig.colorbar(p, cax=cbaxes, orientation='horizontal')

        cbaxes.set_xticklabels([max(0.01, np.round(x, 2)) for x in cbaxes.get_xticks()*max(den)], fontsize=13)
        cbaxes.set_zorder(1)
        axs[z].text(8.2, 5.7, r'$z = {}$'.format(z), fontsize = 18)
        if z != 0:
            axs[z].legend().set_visible(False)

else:

    print ('Not an applicable choice, retry....')
    plt.close()
    sys.exit()


fig.tight_layout()
fig.subplots_adjust(bottom=0.09, left = 0.08, wspace=0, hspace=0)
if inp not in [4, 5]:
    fig.text(0.03, 0.5, ylab, va='center', rotation='vertical', fontsize=26)
else:
    fig.text(0.01, 0.5, ylab, va='center', rotation='vertical', fontsize=26)

fig.text(0.46, 0.04, xlab, va='center', fontsize=26)
fig.savefig(savename+add+'_full_age.pdf')
plt.close()
print ('End time is ', str(datetime.datetime.now()))
