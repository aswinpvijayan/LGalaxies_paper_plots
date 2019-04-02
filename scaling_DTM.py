
import time
import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
from uncertainties import unumpy
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import get_
import gc
import seaborn as sns
sns.set_context("paper")
sns.set_style(style='white')
from scipy.optimize import curve_fit

"""
    Figure 7 in paper
"""

h = 0.673        
redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])

def get_vals(z):
    
    snap = snapnumMR[np.where(redshift == str(z))[0][0]]
    files = '../Dust_output/MR/SA_output_*'
    Mstar = (get_.get_var(files, 'StellarMass', snap)*1e10)/h
    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    Mdust = Mdust1 + Mdust2
    ok = np.logical_and(Mstar > 10**9.0, Mdust > 0)
    
    Mstar = Mstar[ok]
    Mdust = Mdust[ok]
    Age = get_.get_var(files, 'MassWeightAge', snap)[ok]
    mu = get_.get_var(files, 'mu_gas', snap)[ok]
    Type = get_.get_var(files, 'Type', snap)[ok]
    Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
    Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
    Mcg = Mcg1 + Mcg2
    
    Mmet = np.nansum(Mcg[:,2:], axis = 1)
    Mcg = np.nansum(Mcg, axis = 1)
    
    snap = snapnumMRII[np.where(redshift == str(z))[0][0]]
    files = '../Dust_output/MRII/SA_output_*'
    tmp = (get_.get_var(files, 'StellarMass', snap)*1e10)/h
    Mdust1 = get_.get_var(files, 'DustColdGasDiff_elements', snap)
    Mdust2 = get_.get_var(files, 'DustColdGasClouds_elements', snap)
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    tmp1 = Mdust1 + Mdust2
    ok = np.logical_and(tmp > 10**7.5, np.logical_and(tmp < 10**9.0, tmp1 > 0))
    
    Mstar = np.append(Mstar, tmp[ok])
    Mdust = np.append(Mdust, tmp1[ok])
    Age = np.append(Age, get_.get_var(files, 'MassWeightAge', snap)[ok])
    mu = np.append(mu, get_.get_var(files, 'mu_gas', snap)[ok])
    Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
    Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
    Type = np.append(Type, get_.get_var(files, 'Type', snap)[ok])
    
    Mmet = np.append(Mmet, np.nansum(Mcg1[:,2:] + Mcg2[:,2:], axis = 1))
    Mcg = np.append(Mcg, np.nansum(Mcg1 + Mcg2, axis = 1))
    
    Z = Mmet/Mcg
    DTM = Mdust/(Mmet)
    ok = np.where(Mcg > 1e6)
    
    Mdust = Mdust[ok] 
    DTM = DTM[ok]
    Z = Z[ok]
    Mstar = Mstar[ok]
    Age = Age[ok]
    Type = Type[ok]
    mu = mu[ok]
    
    out = get_.remove_(np.array([Mdust, DTM, Z, Mstar, Age, Type, mu]), np.array([5]))
    out = out[(out[5] == 0)]
    out = np.array(out).T
    
    return out[0], out[1], out[2], out[3], out[4], out[5]

for z in range(0, 9):
    
    if z == 0:
        
        try:    
            tmp = np.load('data/DTM_fit_z{}.npz'.format(z))
            Mdust, DTM, Z, Mstar, Age, mu = tmp['Mdust'], tmp['DTM'], tmp['Z'], tmp['Mstar'], tmp['Age'], tmp['mu']
        except:
            Mdust, DTM, Z, Mstar, Age, mu = get_vals(z)
            np.savez('data/DTM_fit_z{}'.format(z), Mdust = Mdust, DTM = DTM, Z = Z, Mstar = Mstar, Age = Age, mu = mu)
        df = pd.DataFrame({'Mdust': Mdust, 'DTM': DTM, 'Z': Z, 'z': np.ones(len(Z))*z, 'Mstar': Mstar, 'Age': Age, 'mu': mu})
    else:
        try:    
            tmp = np.load('data/DTM_fit_z{}.npz'.format(z))
            Mdust, DTM, Z, Mstar, Age, mu = tmp['Mdust'], tmp['DTM'], tmp['Z'], tmp['Mstar'], tmp['Age'], tmp['mu']
        except:
            Mdust, DTM, Z, Mstar, Age, mu = get_vals(z)
            np.savez('data/DTM_fit_z{}'.format(z), Mdust = Mdust, DTM = DTM, Z = Z, Mstar = Mstar, Age = Age, mu = mu)
        data = pd.DataFrame({'Mdust': Mdust, 'DTM': DTM, 'Z': Z, 'z': np.ones(len(Z))*z, 'Mstar': Mstar, 'Age': Age, 'mu': mu})
        df = df.append(data, ignore_index = True)

print ('Data collection complete')

def model(z, theta):
    
    x, y = z
    a, b, c, d, e = theta
    Zsun = 0.0134
    tau = 5e-5/((10**a)*x*Zsun)
    return a + np.log10(1. + b*(np.exp(-c*(x**d)*((y/tau)**e))))

def model_cf(z, a, b, c, d, e):

    Zsun = 0.0134
    tau = 5e-5/((10**a)*z[0]*Zsun)
    return a + np.log10(1. + b*(np.exp(-c*(z[0]**d)*((z[1]/tau)**e))))
    
Zsun = 0.0134
Mdust = np.array(df['Mdust'])
mu = np.array(df['mu'])
DTM = np.array(df['DTM'])
Z = np.array(df['Z'])/Zsun
z = np.array(df['z'])
Age = np.array(df['Age'])
Mstar = np.array(df['Mstar'])

x_obs = [Z, Age]
y_obs = np.log10(DTM)

p0 = [-1.93735568, 30.88799955,  0.60868264,  0.83118602, -1.37868402]

popt, pcov = curve_fit(model_cf, x_obs, y_obs, p0, method = 'lm', sigma = np.ones(len(DTM))*0.1)
print (popt)
print (np.sqrt(pcov.diagonal()))

sol = unumpy.uarray(popt, np.sqrt(pcov.diagonal()))
D0 = 10**(sol[0])
D1 = (sol[1] + 1.)*(10**(sol[0]))
alpha = sol[2]
beta = sol[3]
gamma = sol[4]

print ('########   Output   ########')
print ('D0 = {}'.format(10**(sol[0])))
print ('D1 = {}'.format((sol[1] + 1.)*(10**(sol[0]))))
print ('alpha = {}'.format(sol[2]))
print ('beta = {}'.format(sol[3]))
print ('gamma = {}'.format(sol[4]))


Zsun = 0.0134
Mdust = np.array(df['Mdust'])
#Mcg = np.array(df['Mcg'])
DTM = np.array(df['DTM'])
Z = np.array(df['Z'])/Zsun
z = np.array(df['z'])
Mstar = np.array(df['Mstar'])
Age = np.array(df['Age'])

median_params = popt#[ -1.93741322  30.33397877   0.59917642   0.75647007  -1.36296454]

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

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(6, 6), sharex=True, sharey=True, facecolor='w', edgecolor='k')

points = np.linspace(min(np.log10(DTM)), max(np.log10(DTM))+2, 20)

axs.plot(points, points, ls = 'dashed', color = 'black', lw = 2)
    
ok = np.where(z < 3)
x, y = np.log10(DTM[ok]), model([Z[ok], Age[ok]], median_params)
p = get_contours(x, y, axs, 'blue')

ok = np.logical_and(z >= 3, z < 6)
x, y = np.log10(DTM[ok]), model([Z[ok], Age[ok]], median_params)
get_contours(x, y, axs, 'green')

ok = np.logical_and(z >= 6, z < 9)
x, y = np.log10(DTM[ok]), model([Z[ok], Age[ok]], median_params)
get_contours(x, y, axs, 'red')

axs.set_ylim((-2.4, 0.1))
axs.set_xlim((-2.4, 0.1))
for label in (axs.get_xticklabels() + axs.get_yticklabels()):
    label.set_fontsize(15)
axs.grid()

axs.text(-2.2, -0.10, r'z $\leq$ 2', color = 'blue', fontsize = 15)
axs.text(-2.3, -0.25, r'3 $\leq$ z $\leq$ 5', color = 'green', fontsize = 15)
axs.text(-2.3, -0.40, r'6 $\leq$ z $\leq$ 8', color = 'red', fontsize = 15)

fig.tight_layout()
fig.subplots_adjust(bottom=0.13, left = 0.15, wspace=0, hspace=0)
axs.set_xlabel(r'$\mathrm{log}_{10}(\mathrm{DTM})$', fontsize = 18)
axs.set_ylabel(r'$\mathrm{Fit(Z,}$ $\mathrm{Age)}$', fontsize = 18)

fig.savefig('plot_fit_allontop.pdf')    

plt.close()
