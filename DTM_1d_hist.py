import time
import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import get_
import seaborn as sns
sns.set_context("paper")
sns.set_style(style='white')

"""
    Figure 4 in paper
"""

h = 0.673        
redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])

def get_DTM(z):

    files = '../Dust_output/MR/SA_output_*'
    snap = snapnumMR[np.where(redshift == str(z))[0][0]]

    Mstar = (get_.get_var(files, 'StellarMass', snap)*1e10)/0.673
    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    Mdust = Mdust1 + Mdust2
    ok = np.logical_and(Mstar > 10**9.0, Mdust > 0)

    Mstar = Mstar[ok]
    Mdust = Mdust[ok]
    Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
    Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
    Mcg = Mcg1 + Mcg2
    
    Mmet = np.nansum(Mcg[:,2:], axis = 1)

    snap = snapnumMRII[np.where(redshift == str(z))[0][0]]
    files = '../Dust_output/MRII/SA_output_*'
    tmp = (get_.get_var(files, 'StellarMass', snap)*1e10)/0.673
    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    tmp1 = Mdust1 + Mdust2
    ok = np.logical_and(tmp > 10**7.0, np.logical_and(tmp < 10**9.0, tmp1 > 0))

    Mstar = np.append(Mstar, tmp[ok])
    Mdust = np.append(Mdust, tmp1[ok])
    Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok] + get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]

    Mmet = np.append(Mmet, np.nansum(Mcg1[:,2:], axis = 1))
    Mcg = np.append(Mcg, Mcg1, axis = 0)

    DTM = Mdust/(Mmet)
    ok = np.where(Mcg > 1e6)[0]
    return DTM[ok]

fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize=(8, 8), sharex=True, sharey=False, facecolor='w', edgecolor='k')
axs = axs.ravel()
bins = np.linspace(-3, 0, 30)
for z in range(0, 9):
    
    try:
        tmp = np.load('data/DTM_fit_z{}.npz'.format(z))
        Mdust, DTM, Z, Mstar, Age = tmp['Mdust'], tmp['DTM'], tmp['Z'], tmp['Mstar'], tmp['Age']
    except:
        DTM = np.log10(get_DTM(z))
    
    DTM = np.log10(DTM)
    weights = np.ones_like(DTM)/float(len(DTM))
    hist, edges = np.histogram(DTM, bins = bins)
    frac = hist/(len(DTM)*(edges[1]-edges[0]))
    left,right = edges[:-1],edges[1:]
    X = np.array([left,right]).T.flatten()
    Y = np.array([frac,frac]).T.flatten()
    axs[z].plot(X, np.log10(Y))
    axs[z].set_xlim([-2.6,0.1])
    axs[z].set_xticks((-2, -1, 0))
    if z in range(3,9):
        axs[z].set_ylim(np.log10([0.05, 5]))
        axs[z].set_yticks(np.log10([0.1, 0.3, 0.7, 1.5, 3]))
        axs[z].set_yticklabels([])
    else:
        axs[z].set_ylim(np.log10([0.1, 10]))
        axs[z].set_yticks(np.log10([0.5, 1, 3, 6]))
        axs[z].set_yticklabels([])
        
    if z in [3, 6]:
        axs[z].set_yticklabels([0.1, 0.3, 0.7, 1.5, 3])
    
    axs[z].grid()
    axs[z].text(-2.5, np.log10(2), r'$z = {}$'.format(z), fontsize = 18)
    for label in (axs[z].get_xticklabels() + axs[z].get_yticklabels()):
        label.set_fontsize(13)

axs[0].set_yticklabels([0.5, 1, 3, 6])    
    
fig.subplots_adjust(wspace=0, hspace=0) 
axs[3].set_ylabel('fraction/dex', fontsize = 16)
axs[7].set_xlabel(r'$\mathrm{log}_{10}(\mathrm{DTM})$', fontsize = 18)
plt.savefig('DTM_fracs.pdf')    
plt.close()

