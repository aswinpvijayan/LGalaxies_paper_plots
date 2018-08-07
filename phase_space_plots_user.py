#!/usr/bin/env python

import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import get_
import mbb 
import gc
import seaborn as sns
sns.set_context("paper")


        ###################################################################################
        
redshift = np.array(['0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8',  '9', '10', '11', '12', '13', '14'])
snapnumMR =  np.array(['58', '38', '30', '25', '22', '19', '17', '15', '13', '12', '11', '10', '9', '8', '7'])
snapnumMRII = np.array(['62', '42', '34', '29', '26', '23', '21', '19', '17', '16', '15', '14',	'13', '12', '11', '10'])

snaps = np.array([snapnumMR, snapnumMRII])
sims = np.array(['MR', 'MRII'])

def outputs(x, y, sSFR, Type, z):
    
    """
    #Function remove_ is for removing any nan's, inifnities and other non-physical values, can specify if the zero values in the array should not be removed, by giving the array number as argument. It then selects the central halos which are above a particular sSFR. Returns the number density in a pixel and the median (in 13 bins in the x axis in logarithmic space) with the values of the 16th percentile and the 84th percentile.
    """
    
    out = get_.remove_(np.array([x, y, sSFR, Type]), np.array([3]))  
    #out = out[(out[2] > get_.sSFR_cut(z)) & (out[3] == 0)]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    this_sSFR = out[2]
    this_Type = out[3]
    
    del out, x, y, sSFR, Type
    
    
    xx, yy, yy_up, yy_low = get_.get_median(thisx, thisy, n = 25)
    
    den = get_.get_density(thisx, thisy)
    
    return thisx, thisy, den, xx, yy, yy_up, yy_low


def plot_Mstar_Mdust(files, z, axs, snapnum, i):
    
    from obs_plots import DM_obs 
    
    if i in [0, 1]:
    
        add = sims[i] 
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        Type = get_.get_var(files[i], 'Type', snap)
        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
        Mdust1 = np.nansum(Mdust1, axis = 1)  
        Mdust2 = np.nansum(Mdust2, axis = 1) 
        Mdust = Mdust1 + Mdust2
        SFR = get_.get_var(files[i], 'Sfr', snap)
    
    else:
        
        add = 'MR_MRII'
        for j in [0, 1]:
            
            snap = snapnum[j][np.where(redshift == str(z))[0][0]]
            if j == 0:
                
                Mstar = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(Mstar >= 10**8.9)
                Mstar = Mstar[ok]
                Type = get_.get_var(files[j], 'Type', snap)[ok]
                Mdust1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                Mdust2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                Mdust1 = np.nansum(Mdust1, axis = 1)  
                Mdust2 = np.nansum(Mdust2, axis = 1) 
                Mdust = Mdust1 + Mdust2
                SFR = get_.get_var(files[j], 'Sfr', snap)[ok]
                
            else:
                
                tmp = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(tmp < 10**8.9)
                Mstar = np.append(Mstar, tmp[ok])
                Type = np.append(Type, get_.get_var(files[j], 'Type', snap)[ok])
                tmp1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                tmp1 = np.nansum(tmp1, axis = 1)  
                tmp2 = np.nansum(tmp2, axis = 1) 
                Mdust = np.append(Mdust, tmp1 + tmp2)
                SFR = np.append(SFR, get_.get_var(files[j], 'Sfr', snap)[ok])
    
    sSFR = SFR/Mstar
    
    x, y, den, xx, yy, yy_up, yy_low = outputs(Mstar, Mdust, sSFR, Type, z)
    
    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    DM_obs(axs, z)   #Plotting the observational data points
    
    plot_figure(axs, z, x, y, den, xx, yy, yy_up, yy_low)
    
    xlim = [7.5,11.7]
    ylim = [2.5,10.5]
    xticks = [7, 8, 9, 10, 11]
    axs.text(8.65, 9.5, r'$z = {}$'.format(z), fontsize = 18)
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    return add

def plot_O_H_vs_Dust(files, z, axs, snapnum, i):

    from obs_plots import D_Met_obs
    
    if i in [0, 1]:
        
        add = sims[i]
        
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Type = get_.get_var(files[i], 'Type', snap)
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)
        
        met = (Mcg1[:,4]+Mcg2[:,4])/(15.9994*(Mcg1[:,0]+Mcg2[:,0]))  #O/H ratio
        
        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
        Mdust1 = np.nansum(Mdust1, axis = 1)  
        Mdust2 = np.nansum(Mdust2, axis = 1) 
        Mdust = Mdust1 + Mdust2    
        
        SFR = get_.get_var(files[i], 'Sfr', snap)
    
    else:
        
        add = 'MR_MRII'
        for j in [0, 1]:
            
            snap = snapnum[j][np.where(redshift == str(z))[0][0]]
            if j == 0:
                
                Mstar = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(Mstar >= 10**8.9)
                
                Mstar = Mstar[ok]
                Type = get_.get_var(files[j], 'Type', snap)[ok]
                
                Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
                Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
                met = (Mcg1[:,4]+Mcg2[:,4])/(15.9994*(Mcg1[:,0]+Mcg2[:,0]))
                
                Mdust1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                Mdust2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                Mdust1 = np.nansum(Mdust1, axis = 1)  
                Mdust2 = np.nansum(Mdust2, axis = 1) 
                Mdust = Mdust1 + Mdust2
                SFR = get_.get_var(files[j], 'Sfr', snap)[ok]
            
            else:
            
                tmp = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(tmp < 10**8.9)
                Mstar = np.append(Mstar, tmp[ok])
                Type = np.append(Type, get_.get_var(files[j], 'Type', snap)[ok])
                
                tmp1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
                tmp = (tmp1[:,4]+tmp2[:,4])/(15.9994*(tmp1[:,0]+tmp2[:,0]))
                met = np.append(met, tmp)
                
                tmp1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                tmp1 = np.nansum(tmp1, axis = 1)  
                tmp2 = np.nansum(tmp2, axis = 1) 
                Mdust = np.append(Mdust, tmp1 + tmp2)
                SFR = np.append(SFR, get_.get_var(files[j], 'Sfr', snap)[ok])
    
    sSFR = SFR/Mstar
    
    x, y, den, xx, yy, yy_up, yy_low = outputs(met, Mdust, sSFR, Type, z)
    
    x = 12 + np.log10(x)
    y = np.log10(y)
    xx = 12 + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    D_Met_obs(axs, z)    #Plotting the observational data points
    plot_figure(axs, z, x, y, den, xx, yy, yy_up, yy_low)
    
    xlim = [7,10.5]
    ylim = [2.5,10.5]
    xticks = [7, 8, 9, 10]
    axs.text(7.75, 9, r'$z = {}$'.format(z), fontsize = 18)
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    return add
    
def plot_DGR_Mstar(files, z, axs, snapnum, i):

    from obs_plots import DG_Mstar_obs
    
    if i in [0, 1]:
        
        add = sims[i]
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Type = get_.get_var(files[i], 'Type', snap)
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)
        Mcg = Mcg1 + Mcg2
        Mcg = np.sum(Mcg, axis = 1) #Total cold gas mass
        
        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
        Mdust = Mdust1 + Mdust2    
        Mdust = np.sum(Mdust, axis = 1)  #Total dust mass
        
        Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
        
        SFR = get_.get_var(files[i], 'Sfr', snap)
        
    else:
    
        add = 'MR_MRII'
        for j in [0, 1]:
            
            snap = snapnum[j][np.where(redshift == str(z))[0][0]]
            if j == 0:
                
                Mstar = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(Mstar >= 10**8.9)
                
                Mstar = Mstar[ok]
                Type = get_.get_var(files[j], 'Type', snap)[ok]
                
                Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
                Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
                Mcg = np.nansum(Mcg1, axis = 1) + np.nansum(Mcg2, axis = 1)
                
                Mdust1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                Mdust2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                Mdust1 = np.nansum(Mdust1, axis = 1)  
                Mdust2 = np.nansum(Mdust2, axis = 1) 
                Mdust = Mdust1 + Mdust2
                
                Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
                
                SFR = get_.get_var(files[j], 'Sfr', snap)[ok]
            
            else:
            
                tmp = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(tmp < 10**8.9)
                Mstar = np.append(Mstar, tmp[ok])
                Type = np.append(Type, get_.get_var(files[j], 'Type', snap)[ok])
                
                tmp1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
                tmp = np.nansum(tmp1, axis = 1) + np.nansum(tmp2, axis = 1)
                
                tmp1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                tmp1 = np.nansum(tmp1, axis = 1)  
                tmp2 = np.nansum(tmp2, axis = 1) 
                
                Mratio = np.append(Mratio, (tmp1+tmp2)/tmp)
                
                SFR = np.append(SFR, get_.get_var(files[j], 'Sfr', snap)[ok])
                
    sSFR = SFR/Mstar
    
    x, y, den, xx, yy, yy_up, yy_low = outputs(Mstar, Mratio, sSFR, Type, z)
    
    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    DG_Mstar_obs(axs, z) #Plotting the observational data points
    plot_figure(axs, z, x, y, den, xx, yy, yy_up, yy_low)
    
    xlim = [8,11.5]
    ylim = [-7,-0.3]
    xticks = [8, 9, 10, 11]
    axs.text(8.75, -1, r'$z = {}$'.format(z), fontsize = 18)
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    return add
    
def plot_O_H_vs_DGR(files, z, axs, snapnum, i):

    from obs_plots import DG_met_obs
    
    if i in [0, 1]:
        
        add = sims[i]
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Type = get_.get_var(files[i], 'Type', snap)
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)
        Mcg = Mcg1 + Mcg2
        
        met = (Mcg[:,4])/(15.9994*(Mcg[:,0]))  #O/H ratio
        
        Mcg = np.sum(Mcg, axis = 1) #Total cold gas mass
        
        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
        Mdust1 = np.nansum(Mdust1, axis = 1)  
        Mdust2 = np.nansum(Mdust2, axis = 1) 
        Mdust = Mdust1 + Mdust2
        
        Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
        
        SFR = get_.get_var(files[i], 'Sfr', snap)
    
    else:
        
        add = 'MR_MRII'
        for j in [0, 1]:
            
            snap = snapnum[j][np.where(redshift == str(z))[0][0]]
            if j == 0:
                
                Mstar = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(Mstar >= 10**8.9)
                
                Mstar = Mstar[ok]
                Type = get_.get_var(files[j], 'Type', snap)[ok]
                
                Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
                Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
                Mcg = np.nansum(Mcg1, axis = 1) + np.nansum(Mcg2, axis = 1)
                
                Mdust1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                Mdust2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                Mdust1 = np.nansum(Mdust1, axis = 1)  
                Mdust2 = np.nansum(Mdust2, axis = 1) 
                Mdust = Mdust1 + Mdust2
                
                Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
                
                SFR = get_.get_var(files[j], 'Sfr', snap)[ok]
            
            else:
            
                tmp = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(tmp < 10**8.9)
                Mstar = np.append(Mstar, tmp[ok])
                Type = np.append(Type, get_.get_var(files[j], 'Type', snap)[ok])
                
                tmp1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
                tmp = (tmp1[:,4]+tmp2[:,4])/(15.9994*(tmp1[:,0]+tmp2[:,0]))
                met = np.append(met, tmp)
                
                tmp1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                tmp1 = np.nansum(tmp1, axis = 1)  
                tmp2 = np.nansum(tmp2, axis = 1) 
                Mratio = np.append(Mratio, (tmp1+tmp2)/tmp)
                
                SFR = np.append(SFR, get_.get_var(files[j], 'Sfr', snap)[ok])
                
    sSFR = SFR/Mstar
    
    x, y, den, xx, yy, yy_up, yy_low = outputs(met, Mratio, sSFR, Type, z)
        
    x = 12 + np.log10(x)
    y = np.log10(y)
    xx = 12 + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    plot_figure(axs, z, x, y, den, xx, yy, yy_up, yy_low)
    DG_met_obs(axs, z)   #Plotting the observational data points
    
    xlim = [5.5,10]
    ylim = [-7,-0.1]
    xticks = [6, 7, 8, 9, 10]
    axs.text(5.75, -1, r'$z = {}$'.format(z), fontsize = 18)
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    
    
def plot_Mstar_DGR(files, z, axs, snapnum, i):
    
    if i in [0, 1]:
        
        add = sims[i]
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Type = get_.get_var(files[i], 'Type', snap)
        
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)
        Mcg = Mcg1[:,0] + Mcg2[:,0]
        #Mcg = np.sum(Mcg, axis = 1) #Total cold gas mass
        
        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
        Mdust = Mdust1 + Mdust2    
        Mdust = np.sum(Mdust, axis = 1)  #Total dust mass
        
        Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
        
        SFR = get_.get_var(files[i], 'Sfr', snap)
    
    else:
        
        add = 'MR_MRII'
        for j in [0, 1]:
            
            snap = snapnum[j][np.where(redshift == str(z))[0][0]]
            if j == 0:
                
                Mstar = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(Mstar >= 10**8.9)
                
                Mstar = Mstar[ok]
                Type = get_.get_var(files[j], 'Type', snap)[ok]
                
                Mcg1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
                Mcg2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
                met = (Mcg1[:,4]+Mcg2[:,4])/(15.9994*(Mcg1[:,0]+Mcg2[:,0]))
                
                Mdust1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                Mdust2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                Mdust1 = np.nansum(Mdust1, axis = 1)  
                Mdust2 = np.nansum(Mdust2, axis = 1) 
                Mdust = Mdust1 + Mdust2
                
                Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
                
                SFR = get_.get_var(files[j], 'Sfr', snap)[ok]
            
            else:
            
                tmp = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(tmp < 10**8.9)
                Mstar = np.append(Mstar, tmp[ok])
                Type = np.append(Type, get_.get_var(files[j], 'Type', snap)[ok])
                
                tmp1 = get_.get_var(files, 'ColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files, 'ColdGasClouds_elements', snap)[ok]
                tmp = np.nansum(tmp1, axis = 1) + np.nansum(tmp2, axis = 1)
                
                tmp1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                tmp1 = np.nansum(tmp1, axis = 1)  
                tmp2 = np.nansum(tmp2, axis = 1) 
                Mratio = np.append(Mratio, (tmp1+tmp2)/tmp)
                
                SFR = np.append(SFR, get_.get_var(files[j], 'Sfr', snap)[ok])
    
    
    sSFR = SFR/Mstar
    
    x, y, den, xx, yy, yy_up, yy_low = outputs(Mstar, Mratio, sSFR, Type, z)
        
    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    plot_figure(axs, z, x, y, den, xx, yy, yy_up, yy_low)
    
    xlim = [6.5,11.5]
    ylim = [-5,-0.2]
    xticks = [7, 8, 9, 10, 11]
    axs.text(8.75, -0.9, r'$z = {}$'.format(z), fontsize = 18)
    
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    return add

def plot_sSFR_MdustMstar(files, z, axs, snapnum, i):
    
    if i in [0, 1]:
        
        add = sims[i]
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Type = get_.get_var(files[i], 'Type', snap)
        
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        
        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
        Mdust = Mdust1 + Mdust2    
        Mdust = np.nansum(Mdust, axis = 1)  #Total dust mass
        
        SFR = get_.get_var(files[i], 'Sfr', snap)
    
    else:
        for j in [0, 1]:
            
            add = 'MR_MRII'
            snap = snapnum[j][np.where(redshift == str(z))[0][0]]
            if j == 0:
                
                Mstar = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(Mstar >= 10**8.9)
                Mstar = Mstar[ok]
                Type = get_.get_var(files[j], 'Type', snap)[ok]
                Mdust1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                Mdust2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                Mdust1 = np.nansum(Mdust1, axis = 1)  
                Mdust2 = np.nansum(Mdust2, axis = 1) 
                Mdust = Mdust1 + Mdust2
                SFR = get_.get_var(files[j], 'Sfr', snap)[ok]
                
            else:
                
                tmp = (get_.get_var(files[j], 'StellarMass', snap)*1e10)/0.673
                ok = np.where(tmp < 10**8.9)
                Mstar = np.append(Mstar, tmp[ok])
                Type = np.append(Type, get_.get_var(files[j], 'Type', snap)[ok])
                tmp1 = get_.get_var(files[j], 'DustColdGasDiff_elements', snap)[ok]
                tmp2 = get_.get_var(files[j], 'DustColdGasClouds_elements', snap)[ok]
                tmp1 = np.nansum(tmp1, axis = 1)  
                tmp2 = np.nansum(tmp2, axis = 1) 
                Mdust = np.append(Mdust, tmp1 + tmp2)
                SFR = np.append(SFR, get_.get_var(files[j], 'Sfr', snap)[ok])
    
    sSFR = SFR/Mstar
    
    Lir =  57.908942450477525*4*np.pi*Mdust   #Conversion factor for a dust temperature of 50K
    
    x, y, den, xx, yy, yy_up, yy_low = outputs(Lir, SFR, sSFR, Type, z)
    
    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    
    plot_figure(axs, z, x, y, den, xx, yy, yy_up, yy_low)
    
    tmp = np.linspace(7.9, 13, 50)
    axs.plot(tmp, -9.6 + tmp - (2./(tmp-7.)), ls = 'dashed', color = 'red') #Median line from Santini et al. 2014
    
    xlim = [8,12.3]
    ylim = [-3,3]
    xticks = [8, 9, 10, 11, 12]
    axs.text(8.5, 2, r'$z = {}$'.format(z), fontsize = 18)
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    return add

def plot_figure(axs, z, x, y, den, xx, yy, yy_up, yy_low):
    
    axs.scatter(x, y, marker = 'o', s=10, c = den, edgecolors='None', alpha = 0.09, cmap = plt.cm.get_cmap('gist_gray'))
    axs.hexbin(x, y, gridsize=500, cmap = plt.cm.get_cmap('gist_gray'), mincnt=10)
    axs.plot(xx, yy, lw = 2, color = 'black')
    axs.plot(xx, yy_up, lw = 2, ls = 'dashed', color = 'grey')
    axs.plot(xx, yy_low, lw = 2, ls = 'dashed', color = 'grey')
    
    del x, y, xx, yy, yy_up, yy_low
    
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(18)
    
    axs.grid(True)
    lgd = axs.legend(fontsize = 15, markerscale=2, loc = 4, numpoints=1, handletextpad=0.005)
    if np.isscalar(lgd):
        lgd.set_zorder(100)


        ####################################################################################

"""
    User inputs for producing the preferred plots: 
    0. Dust mass vs stellar mass plot for z = 0-9
    1. Dust mass vs Metallciity plot for z = 0-9
    2. Dust-to-gas ratio (DGR) vs stellar mass for z = 0-9
    3. DGR vs Metallicity for z = 0-9
    4. Dust-to-metal ratio vs Stellar mass for z = 0-9
    5. SFR vs Lir for z = 0
"""

inp = int(input('Enter the appropriate number for the required phase space plot: \n0 - Dust mass vs stellar mass plot for z = 0-9 \n1 - Dust mass vs Metallciity plot for z = 0-9\n2 - Dust-to-gas ratio (DGR) vs stellar mass for z = 0-9\n3 - DGR vs Metallicity for z = 0-9\n4 - Dust-to-metal ratio vs Stellar mass for z = 0-9\n5 - SFR vs Lir for z = 0\n'))

filesMR = '/lustre/scratch/astro/ap629/Dust_output_8june/MR/SA_output_*'
filesMRII = '/lustre/scratch/astro/ap629/Dust_output_8june/MRII/SA_output_*'
files = np.array([filesMR, filesMRII])

i = int(input('Input 0 to plot MR, 1 for MRII and anything else for plotting both. \nWhen plotting both MR and MRII a stellar mass cut of 10^8.9 is applied\n'))


if inp in [0, 1, 2, 3, 4]:
    
    fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize=(15, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()
else:
    
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(15, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')



if inp == 0:
    
    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(M_{Dust}/(M_{*}))$'
    savename = 'Dust_Stellar_'
    
    for z in range(1, 9):
        add = plot_Mstar_Mdust(files, z, axs[z], snaps, i)
        gc.collect()
        
elif inp == 1:

    xlab = r'$12+log_{10}(O/H)$'
    ylab = r'$log_{10}(M_{Dust}/(M_{\odot}))$'
    savename = 'Dust_metal_'
    for z in range(1, 9):
        add = plot_O_H_vs_Dust(files, z, axs[z], snaps, i)
        gc.collect()
    
elif inp == 2:
    
    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(M_{Dust}/(M_{cold\mathrm{ }gas}))$'
    savename = 'Dust_gas_ratio_'
    for z in range(1, 9):
        add = plot_DGR_Mstar(files, z, axs[z], snaps, i)
        gc.collect()
        
elif inp == 3:

    xlab = r'$12+log_{10}(O/H)$'
    ylab = r'$log_{10}(M_{Dust}/(M_{cold\mathrm{ }gas}))$'
    savename = 'Met_Dust_gas_ratio_'
    for z in range(1, 9):
        add = plot_O_H_vs_DGR(files, z, axs[z], snaps, i)
        gc.collect()
        
elif inp == 4:

    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(M_{Dust}/(M_{Metal}))$'
    savename = 'Dust_metal_ratio_'
    for z in range(1, 9):
        add = plot_Mstar_DGR(files, z, axs[z], snaps, i)
        gc.collect()
        
elif inp == 5:

    xlab = r'$log_{10}(Lir/(L_{\odot}))$'
    ylab = r'$log_{10}(SFR/(M_{\odot}yr^{-1}))$'
    savename = 'Dust_SFR_Lir_'
    add = plot_sSFR_MdustMstar(files, 0, axs, snaps, i)
    gc.collect()
    
else:

    print ('Not an applicable choice, retry....')
    plt.close()
    sys.exit() 


fig.tight_layout()    
fig.subplots_adjust(bottom=0.09, left = 0.08, wspace=0, hspace=0)
fig.text(0.03, 0.5, ylab, va='center', rotation='vertical', fontsize=22)
fig.text(0.5, 0.04, xlab, va='center', fontsize=22)
plt.savefig(savename+add+'.png')
plt.show()

