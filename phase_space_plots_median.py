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
sns.set_style(style='white')
from astropy.cosmology import Planck13

        #####################################################################################
"""
    User inputs for producing the preferred plots:
    0. Dust mass vs stellar mass plot    
    2. Dust-to-gas ratio (DGR) vs stellar mass 
    3. DGR vs Metallicity
    4. Dust-to-metal (DTM) ratio vs Stellar mass   
    5. DTM ratio vs Metallicity
    6. Accretion timescale vs stellar mass plot
    7. Accretion timescale vs metallicity plot

    inp = user input for preferred plot

    z_inp = Select redshift (z) range to plot:
    0. [0, 8]
    1. [9, 13]
    2. [0, 13]

    i = 0 to plot just MR, 1 for MRII and any higher number for plotting both MR and MRII

"""

filesMR = '../Dust_output/MR/SA_output_*'
filesMRII = '../Dust_output/MRII/SA_output_*'
files = np.array([filesMR, filesMRII])

inp, z_inp, i = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])

        ###################################################################################

h = 0.673
redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

norm_z = np.arange(0, 9, 1)
high_z = np.arange(9, 13, 1)
all_z = np.arange(0, 13, 1)

snaps = np.array([snapnumMR, snapnumMRII])
sims = np.array(['MR', 'MRII'])

MR_vol = (480.279/h)**3  #Millennium
MRII_vol = (96.0558/h)**3 #Millennium II

if z_inp == 0:
    z_ = norm_z
    name_z = 'z0_8'
    
elif z_inp == 1:
    z_ = high_z
    name_z = 'z9_12'
    
elif z_inp == 2:
    z_ = all_z
    name_z = 'z0_12'

else:
    print ('Not an applicable redshift choice, retry....')
    sys.exit()

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(6, 6), sharex=True, sharey=True, facecolor='w', edgecolor='k')


if inp == 0:

    xlab = r'$\mathrm{log}_{10}(\mathrm{M}_{*}/\mathrm{M}_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(\mathrm{M}_{\mathrm{Dust}}/M_{*})$'
    savename = 'Dust_Stellar_'

    from Mstar_Mdust import plot_Mstar_Mdust_median
    if z_inp == 0:
        xlim = [7.5,11.9]
        ylim = [1.5,10.8]
        xticks = [8, 9, 10, 11]
    else:
        xlim = [7,11.7]
        ylim = [1.5,6.5]
        xticks = [7, 8, 9, 10, 11]
        
    for z in z_:
        if i <= 1:
            add = plot_Mstar_Mdust_median(files, z, axs, snaps, i, False)
        else:
            add1 = plot_Mstar_Mdust_median(files, z, axs, snaps, 0, True)
            #print ('MR: z={}'.format(z))
            add2 = plot_Mstar_Mdust_median(files, z, axs, snaps, 1, True)
            #print ('MRII: z={}'.format(z))
            add = add1+'_'+add2
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)

elif inp == 1:

    xlab = r'$12+\mathrm{log}_{10}(\mathrm{O/H})$'
    ylab = r'$\mathrm{log}_{10}(\mathrm{M}_{\mathrm{Dust}}/\mathrm{M}_{\odot})$'
    savename = 'Dust_metal_'

    from Met_Mdust import plot_O_H_vs_Dust_median

    for z in z_:
        if i <= 1:
            add = plot_O_H_vs_Dust_median(files, z, axs, snaps, i, False)
        else:
            add = plot_O_H_vs_Dust_median(files, z, axs, snaps, i, True)

elif inp == 2:

    xlab = r'$\mathrm{log}_{10}(\mathrm{M}_{*}/\mathrm{M}_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(\mathrm{M}_{\mathrm{Dust}}/\mathrm{M}_{\mathrm{Cold gas}})$'
    savename = 'DGR_stell_'

    from DGR_Mstar import plot_DGR_Mstar_median

    for z in z_:
        if i <= 1:
            add = plot_DGR_Mstar_median(files, z, axs, snaps, i, False)
        else:
            add1 = plot_DGR_Mstar_median(files, z, axs, snaps, 0, True)
            add2 = plot_DGR_Mstar_median(files, z, axs, snaps, 1, True)
            add = add1+'_'+add2

elif inp == 3:

    xlab = r'$12+\mathrm{log}_{10}(O/H)$'
    ylab = r'$\mathrm{log}_{10}(\mathrm{M}_{\mathrm{Dust}}/\mathrm{M}_{\mathrm{Cold gas}})$'
    savename = 'DGR_met_ratio_'

    from DGR_Met import plot_O_H_vs_DGR_median

    for z in z_:
        if i <= 1:
            add = plot_O_H_vs_DGR_median(files, z, axs, snaps, i, False)
        else:
            add = plot_O_H_vs_DGR_median(files, z, axs, snaps, i, True)

elif inp == 4:

    xlab = r'$\mathrm{log}_{10}(\mathrm{M}_{*}/\mathrm{M}_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(\mathrm{M}_{\mathrm{Dust}}/\mathrm{M}_{\mathrm{Metal}})$'
    savename = 'DTM_stellar_'
    
    if z_inp == 0:
        xlim = [7.5, 11.7]
        ylim = [-2, 0]
        xticks = [8, 9, 10, 11]
    else:
        xlim = [7,11.3]
        ylim = [-2.5, -1.5]
        xticks = [7, 8, 9, 10, 11]

    from DTM_Mstar import plot_Mstar_DTM_median

    for z in z_:
        if i <= 1:
            add = plot_Mstar_DTM_median(files, z, axs, snaps, i, False)
            print ('{}: z={}'.format(sims[i], z))
        else:
            add1 = plot_Mstar_DTM_median(files, z, axs, snaps, 0, True)
            print ('MR: z={}'.format(z))
            add2 = plot_Mstar_DTM_median(files, z, axs, snaps, 1, True)
            print ('MRII: z={}'.format(z))
            add = add1+'_'+add2

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    #axs.set_yticks(yticks)
    
elif inp == 5:

    xlab = r'$12+\mathrm{log}_{10}(O/H)$'
    ylab = r'$\mathrm{log}_{10}(\mathrm{M}_{\mathrm{Dust}}/\mathrm{M}_{\mathrm{Metal}})$'
    savename = 'DTM_met_'

    from DTM_Met import plot_O_H_vs_DTM_median

    for z in z_:
        if i <= 1:
            add = plot_O_H_vs_DTM_median(files, z, axs, snaps, i, False)
            print ('{}: z={}'.format(sims[i], z))
        else:
            add = plot_O_H_vs_DTM_median(files, z, axs, snaps, i, True)
            print ('z={}'.format(z))

elif inp == 6:

    xlab = r'$\mathrm{log}_{10}(M_{*}/M_{\odot})$'
    ylab = r'$\mathrm{log}_{10}(\tau_{\mathrm{acc}})$'
    savename = 'Mstar_tacc_'

    from tacc_Mstar import plot_tacc_Mstar_median

    for j,z in enumerate(z_):

        if i <= 1:
            add = plot_tacc_Mstar_median(files, z, axs, snaps, i, False)
            print ('{}: z={}'.format(sims[i], z))
        else:
            add1 = plot_tacc_Mstar_median(files, z, axs, snaps, 0, True)
            print ('MR: z={}'.format(sims[i], z))
            add2 = plot_tacc_Mstar_median(files, z, axs, snaps, 1, True)
            print ('MRII: z={}'.format(sims[i], z))
            add = add1+'_'+add2

elif inp == 7:

    xlab = r'$12+\mathrm{log}_{10}(\mathrm{O/H})$'
    ylab = r'$\mathrm{log}_{10}(\tau_{\mathrm{acc}})$'
    savename = 'Met_tacc_'

    from tacc_Met import plot_tacc_Met_median

    for j,z in enumerate(z_):

        if i <= 1:
            add = plot_tacc_Met_median(files, z, axs, snaps, i, False)
            print ('{}: z={}'.format(sims[i], z))
        else:
            add = plot_tacc_Met_median(files, z, axs, snaps, i, True)
            print ('z={}'.format(sims[i], z))

else:

    print ('Not an applicable choice, retry....')
    plt.close()
    sys.exit()

fig.tight_layout()
fig.subplots_adjust(bottom=0.15, left = 0.15, wspace=0, hspace=0)
fig.text(0.03, 0.5, ylab, va='center', rotation='vertical', fontsize=22)
fig.text(0.45, 0.04, xlab, va='center', fontsize=22)
plt.savefig(name_z+savename+add+'_median.pdf')
plt.close()
