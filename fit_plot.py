
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
import gc
import seaborn as sns
sns.set_context("paper")
sns.set_style(style='white') 

h = 0.673        
redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])

def get_vals(z):
    
    snap = snapnumMR[np.where(redshift == str(z))[0][0]]
    #file_req = '/lustre/scratch/astro/ap629/Dust_output_17jan/MR/hdf5_mod/lgal_z{}_N*.hdf5'.format(z)
    files = '../Dust_output/MR/SA_output_5.h5'
    Mstar = (get_.get_var(files, 'StellarMass', snap)*1e10)/0.673
    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    Mdust = Mdust1 + Mdust2
    ok = np.logical_and(Mstar > 10**8.7, Mdust > 0)
    
    Mstar = Mstar[ok]
    Mdust = Mdust[ok]
    Age = get_.get_var(files, 'MassWeightAge', snap)[ok]
    Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
    Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
    Mcg = Mcg1 + Mcg2
    
    Mmet = np.nansum(Mcg[:,2:], axis = 1)
    Mcg = np.nansum(Mcg, axis = 1)
    
    snap = snapnumMRII[np.where(redshift == str(z))[0][0]]
    files = '../Dust_output/MRII_sub/SA_output_*'
    tmp = (get_.get_var(files, 'StellarMass', snap)*1e10)/0.673
    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    tmp1 = Mdust1 + Mdust2
    ok = np.logical_and(tmp > 10**7.5, np.logical_and(tmp < 10**8.75, tmp1 > 0))
    
    Mstar = np.append(Mstar, tmp[ok])
    Mdust = np.append(Mdust, tmp1[ok])
    Age = np.append(Age, get_.get_var(files, 'MassWeightAge', snap)[ok])
    Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
    Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
    
    Mmet = np.append(Mmet, np.nansum(Mcg1[:,2:] + Mcg2[:,2:], axis = 1))
    Mcg = np.append(Mcg, np.nansum(Mcg1 + Mcg2, axis = 1))
    
    Z = Mmet/Mcg
    DTM = Mdust/(Mmet)
    ok = np.where(Mcg > 1e6)
    return Mdust[ok], DTM[ok], Z[ok], Mstar[ok], Age[ok]

def get_contours(x, y, axs, colour):
    
    n = len(x)
    x_edges = np.arange(min(x)-0.05, max(x)+0.05, 0.03)
    y_edges = np.arange(min(y)-0.05, max(y)+0.05, 0.03)
    hist, xedges, yedges = np.histogram2d(x, y, bins=(x_edges, y_edges))

    xidx = np.digitize(x, x_edges)-1
    yidx = np.digitize(y, y_edges)-1
    num_pts = (hist[xidx, yidx])
    srt = np.sort(num_pts)
    perc = np.array([0.05, 0.32, 0.5])*n
    perc = perc.astype(int)
    levels = srt[perc]
    
    axs.contour(hist.T, extent=[xedges.min(),xedges.max(),yedges.min(),y_edges.max()], levels = levels, colors = colour, linestyles = ['dotted', 'dashed', 'solid'])
    
    return axs


"""
for z in range(0, 9):
    
    if z == 0:
        
        #Mdust, DTM, Z, Mstar, Age = get_vals(z)
        #np.savez('data/DTM_fit_z{}'.format(z), Mdust = Mdust, DTM = DTM, Z = Z, Mstar = Mstar, Age = Age)
        tmp = np.load('data/DTM_fit_z{}.npz'.format(z))
        Mdust, DTM, Z, Mstar, Age = tmp['Mdust'], tmp['DTM'], tmp['Z'], tmp['Mstar'], tmp['Age']
        df = pd.DataFrame({'Mdust': Mdust, 'DTM': DTM, 'Z': Z, 'z': np.ones(len(Z))*z, 'Mstar': Mstar, 'Age': Age})
    else:
        #Mdust, DTM, Z, Mstar, Age = get_vals(z)
        #np.savez('data/DTM_fit_z{}'.format(z), Mdust = Mdust, DTM = DTM, Z = Z, Mstar = Mstar, Age = Age)
        tmp = np.load('data/DTM_fit_z{}.npz'.format(z))
        Mdust, DTM, Z, Mstar, Age = tmp['Mdust'], tmp['DTM'], tmp['Z'], tmp['Mstar'], tmp['Age']
        data = pd.DataFrame({'Mdust': Mdust, 'DTM': DTM, 'Z': Z, 'z': np.ones(len(Z))*z, 'Mstar': Mstar, 'Age': Age})
        df = df.append(data, ignore_index = True)
"""
print ('Data collection complete')

Zsun = 0.0134
Mdust = np.array(df['Mdust'])
#Mcg = np.array(df['Mcg'])
DTM = np.array(df['DTM'])
Z = np.array(df['Z'])/Zsun
z = np.array(df['z'])
Mstar = np.array(df['Mstar'])
Age = np.array(df['Age'])


def model(z, theta):
    
    x, y = z
    a, b, c, d, e = theta
    Zsun = 0.0134
    tau = 5e-5/((10**a)*x*Zsun)
    return a + np.log10(1. + b*(np.exp(-c*(x**d)*((y/tau)**e))))

median_params = [-1.93735568, 30.88799955,  0.60868264,  0.83118602, -1.37868402]
fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(8, 8), sharex=True, sharey=True, facecolor='w', edgecolor='k')
points = np.linspace(min(np.log10(DTM)), max(np.log10(DTM))+2, 20)

axs.plot(points, points, ls = 'dashed', color = 'black', lw = 2)

axs.text(-2.2, -0.1, r'z $\leq$ 2', color = 'blue', fontsize = 15)
axs.text(-2.3, -0.22, r'3 $\leq$ z $\leq$ 5', color = 'green', fontsize = 15)
axs.text(-2.3, -0.34, r'6 $\leq$ z $\leq$ 8', color = 'red', fontsize = 15)
    
ok = np.where(z < 3)
x, y = np.log10(DTM[ok]), model([Z[ok], Age[ok]], median_params)
p = get_contours(x, y, axs, 'blue')

ok = np.logical_and(z >= 3, z < 6)
x, y = np.log10(DTM[ok]), model([Z[ok], Age[ok]], median_params)
get_contours(x, y, axs, 'green')

ok = np.logical_and(z >= 6, z < 9)
x, y = np.log10(DTM[ok]), model([Z[ok], Age[ok]], median_params)
get_contours(x, y, axs, 'red')

axs.set_ylim((-2.5, 0))
axs.set_xlim((-2.5, 0))
for label in (axs.get_xticklabels() + axs.get_yticklabels()):
    label.set_fontsize(15)
axs.grid()

#fig.subplots_adjust(wspace=0, hspace=0, right=0.85) 
#cbar_ax = fig.add_axes([0.9, 0.2, 0.02, 0.55])
#fig.colorbar(p, cax=cbar_ax)   
#cbar_ax.set_yticklabels([np.round(x, 1) for x in cbar_ax.get_yticks()/max(cbar_ax.get_ylim())], fontsize=11)
#cbar_ax.set_ylabel('Normalised 2D density / Redshift', fontsize = 12)
axs.set_xlabel(r'$\mathrm{log}_{10}(\mathrm{M_{dust}} / \mathrm{M_{met}})$', fontsize = 20)
axs.set_ylabel('Fit(Z, Age)', fontsize = 16)
fig.savefig('plot_fit_allontop.pdf')    

plt.close()
