#!/usr/bin/env python

import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
sys.path.append('../')   
import numpy as np
import pandas as pd
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

    if on == False: 
        add = sims[i]
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
        
        Age = get_.get_var(files[i], 'MassWeightAge', snap)
        
        Type = get_.get_var(files[i], 'Type', snap)
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)
        Mcg = Mcg1 + Mcg2
        met = Mcg[:,4]/(15.9994*Mcg[:,0])  #O/H ratio
        
        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
        Mdust = Mdust1 + Mdust2    
        
        met = met - Mdust[:,4]/(15.9994*Mcg[:,0]) 
        
        Mdust = np.nansum(Mdust, axis = 1)  
        
    else:
            i = 0
            snap = snapnum[i][np.where(redshift == str(z))[0][0]]
            Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
            ok = np.where(Mstar >= 10**8.9)[0]
                
            Mstar = Mstar[ok]
            Age = get_.get_var(files[i], 'MassWeightAge', snap)[ok]
            Type = get_.get_var(files[i], 'Type', snap)[ok]
            Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok]
            Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]
            Mcg = Mcg1 + Mcg2
            
            met = Mcg[:,4]/(15.9994*Mcg[:,0])  #O/H ratio
            
            Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
            Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
            Mdust = Mdust1 + Mdust2    
            
            met = met - Mdust[:,4]/(15.9994*Mcg[:,0]) 
            
            Mdust = np.nansum(Mdust, axis = 1)  
            
            
            i = 1
            snap = snapnum[i][np.where(redshift == str(z))[0][0]]
            tmp = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
            ok = np.logical_and(tmp > 10**6.0, tmp < 10**8.9)
            Mstar = np.append(Mstar, tmp[ok])
            
            Type = np.append(Type, get_.get_var(files[i], 'Type', snap)[ok])
            Age = np.append(Age, get_.get_var(files[i], 'MassWeightAge', snap)[ok])
            Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok]
            Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]
            
            Mcg = Mcg1 + Mcg2
            
            Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
            Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
            tmp = Mdust1 + Mdust2
            
            met = np.append(met, (Mcg[:,4]-tmp[:,4])/(15.9994*Mcg[:,0]))
            
            Mdust = np.append(Mdust, np.nansum(tmp, axis = 1)) 
            
            add = 'MR_MRII'
            
    return add, met, Mdust, Type, Age
    

def plot_O_H_vs_Dust_user(files, z, axs, snapnum, i, on):

    add, met, Mdust, Type, Age = get_vals(files, z, axs, snapnum, i, on)
            
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_user(met, Mdust, Type, z)
    
    x = 12 + np.log10(x)
    y = np.log10(y)
    xx = 12 + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    make_fig.fig_user(axs, z, x, y, xx, yy, yy_up, yy_low, den)
    
    xlim = [7.5,10.5]
    ylim = [2.5,10.5]
    xticks = [8, 9, 10]
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    return add
    
def plot_O_H_vs_Dust_median(files, z, axs, snapnum, i, on):

    add, met, Mdust, Type, Age = get_vals(files, z, axs, snapnum, i, on)
            
    xx, yy, yy_up, yy_low = create_out.out_median(met, Mdust, Type, z)

    xx = 12 + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    make_fig.fig_median(axs, z, xx, yy, yy_up, yy_low, 1)
    
    xlim = [6,10.5]
    ylim = [0.1,6.2]
    xticks = [6, 7, 8, 9, 10]
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    return add
    
def plot_O_H_vs_Dust_age(files, z, axs, snapnum, i, on):

    add, met, Mdust, Type, Age = get_vals(files, z, axs, snapnum, i, on)
            
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_age(met, Mdust, Type, Age, z)

    x = 12 + np.log10(x)
    y = np.log10(y)
    xx = 12 + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    p = make_fig.fig_age(axs, z, x, y, xx, yy, yy_up, yy_low, den)
    
    xlim = [7.5,10.5]
    ylim = [2.5,10.5]
    xticks = [8, 9, 10]
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    return add, p, den
