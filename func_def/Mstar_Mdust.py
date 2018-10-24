#!/usr/bin/env python

import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
sys.path.append('../')
import get_
import create_out
import make_fig

h = 0.673        
redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])
sims = np.array(['MR', 'MRII'])

MR_vol = (480.279/h)**3  #Millennium
MRII_vol = (96.0558/h)**3 #Millennium II


def get_vals(files, z, axs, snapnum, i, on):
    
    snap = snapnum[i][np.where(redshift == str(z))[0][0]]
    Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
    if on and i == 0:
        ok = np.where(Mstar >= 10**8.9)[0]
    elif on and i == 1:
        ok = np.logical_and(Mstar > 10**7.5, Mstar < 10**8.9)
    else:
        ok = np.array([True]*len(Mstar))
        
    Mstar = Mstar[ok]
    Type = get_.get_var(files[i], 'Type', snap)[ok]
    Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
    Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    Mdust = Mdust1 + Mdust2
    
    return Mstar, Mdust, Type
    

def plot_Mstar_Mdust_user(files, z, axs, snapnum, i, on):
    
    add = sims[i] 
    
    Mstar, Mdust, Type = get_vals(files, z, axs, snapnum, i, on)
    
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_user(Mstar, Mdust, Type, z)
    
    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    make_fig.fig_user(axs, z, x, y, xx, yy, yy_up, yy_low, den)
    
    xlim = [7.5,11.9]
    ylim = [1.5,10.8]
    xticks = [8, 9, 10, 11]
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    return add
    

def plot_Mstar_Mdust_median(files, z, axs, snapnum, i, on):
    
    Mstar, Mdust, Type = get_vals(files, z, axs, snapnum, i, on)
    
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_median(Mstar, Mdust, Type, z)
    
    xx, yy, yy_up, yy_low = outputs(Mstar, Mdust, Type, z)
    
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    make_fig.fig_med(axs, z, xx, yy, yy_up, yy_low, i)
    
    xlim = [7.5,12.7]
    ylim = [1.5,10.5]
    xticks = [8, 9, 10, 11, 12]
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    return add
    

def plot_Mstar_Mdust_age(files, z, axs, snapnum, i, on):
    
    add = sims[i] 
    snap = snapnum[i][np.where(redshift == str(z))[0][0]]
    
    try:
        data = np.load('data/Mstar_{}_snap_{}.npz'.format(add, snap))
        Mstar = data['Mstar']
        data = np.load('data/Type_{}_snap_{}.npz'.format(add, snap))
        Type = data['Type']
        data = np.load('data/Age_{}_snap_{}.npz'.format(add, snap))
        Age = data['Age']
        data = np.load('data/DustColdGasDiff_elements_{}_snap_{}.npz'.format(add, snap))
        Mdust1 = data['DustColdGasDiff_elements']
        data = np.load('data/DustColdGasClouds_elements_{}_snap_{}.npz'.format(add, snap))
        Mdust2 = data['DustColdGasClouds_elements']
    
    except:
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
        
        if on and i == 0:
            ok = np.where(Mstar >= 10**8.9)[0]
        elif on and i == 1:
            ok = np.logical_and(Mstar > 10**7.5, Mstar < 10**8.9)
        else:
            ok = np.array([True]*len(Mstar))
        
        Mstar = Mstar[ok]
        np.savez_compressed('data/Mstar_{}_snap_{}'.format(add, snap), Mstar=Mstar)
        Type = get_.get_var(files[i], 'Type', snap)[ok]
        np.savez_compressed('data/Type_{}_snap_{}'.format(add, snap), Type=Type)
        Age = get_.get_var(files[i], 'MassWeightAge', snap)[ok]
        np.savez_compressed('data/Age_{}_snap_{}'.format(add, snap), Age=Age)
        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
        np.savez_compressed('data/DustColdGasDiff_elements_{}_snap_{}'.format(add, snap), DustColdGasDiff_elements=Mdust1)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
        np.savez_compressed('data/DustColdGasClouds_elements_{}_snap_{}'.format(add, snap), DustColdGasClouds_elements=Mdust2)
    
    Mdust1 = np.nansum(Mdust1, axis = 1)  
    Mdust2 = np.nansum(Mdust2, axis = 1) 
    Mdust = Mdust1 + Mdust2
    
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_age(Mstar, Mdust, Type, Age, z)
    
    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    p = make_fig.fig_age(axs, z, x, y, xx, yy, yy_up, yy_low, den)
    
    xlim = [7.5,11.9]
    ylim = [1.5,10.8]
    xticks = [8, 9, 10, 11]
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    return add, p, den
