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

    if on:
        i = 0
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
        ok = np.where(Mstar >= 10**9)[0]
        Mstar = Mstar[ok]
        Type = get_.get_var(files[i], 'Type', snap)[ok]
        Age = get_.get_var(files[i], 'MassWeightAge', snap)[ok]
        Mcg = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok] + get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]

        met = (Mcg[:,4])/(15.9994*(Mcg[:,0]))  #O/H ratio

        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
        Mdust = Mdust1 + Mdust2
        met = met - (Mdust[:,4])/(15.9994*(Mcg[:,0]))

        Mcg = np.sum(Mcg, axis = 1) #Total cold gas mass
        Mdust = np.nansum(Mdust, axis = 1)

        Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.

        i = 1
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        tmp = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
        ok = np.logical_and(tmp > 10**7.5, tmp < 10**8.98)
        Mstar = np.append(Mstar, tmp[ok])
        Type = np.append(Type, get_.get_var(files[i], 'Type', snap)[ok])
        Age = np.append(Age, get_.get_var(files[i], 'MassWeightAge', snap)[ok])
        tmp1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok] + get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]

        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)[ok]
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)[ok]
        tmp2 = Mdust1 + Mdust2

        met = np.append(met, (tmp1[:,4]-tmp2[:,4])/(15.9994*(tmp1[:,0])))

        tmp1 = np.nansum(tmp1, axis = 1)
        tmp2 = np.nansum(Mdust1 + Mdust2, axis = 1)

        Mdust = np.append(Mdust, tmp2)

        Mratio = np.append(Mratio, tmp2/tmp1)
        
        Mcg = np.append(Mcg, tmp1)
        ok = np.logical_and(Mcg > 1e6, 12.+np.log10(met) > 6.)  #Only selecting galaxies with cold gas mass above this value, so the galaxies we are looking at are realistic for the dust content usually seen

        add = 'MR_MRII'
    else:
        add = sims[i]
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/h
        
        Type = get_.get_var(files[i], 'Type', snap)
        Age = get_.get_var(files[i], 'MassWeightAge', snap)
        Mcg = get_.get_var(files[i], 'ColdGas_elements', snap)

        met = (Mcg[:,4])/(15.9994*(Mcg[:,0]))  #O/H ratio

        Mcg = np.nansum(Mcg, axis = 1) #Total cold gas mass

        Mdust1 = get_.get_var(files[i], 'DustColdGasDiff_elements', snap)
        Mdust2 = get_.get_var(files[i], 'DustColdGasClouds_elements', snap)
        Mdust1 = np.nansum(Mdust1, axis = 1)
        Mdust2 = np.nansum(Mdust2, axis = 1)
        Mdust = Mdust1 + Mdust2

        Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
        
        ok = np.logical_and(Mcg > 1e6, 12.+np.log10(met) > 6.)  #Only selecting galaxies with cold gas mass above this value, so the galaxies we are looking at are realistic for the dust content usually seen
        
    return add, met[ok], Mratio[ok], Type[ok], Age[ok]


def plot_O_H_vs_DGR_user(files, z, axs, snapnum, i, on):

    add, met, Mratio, Type, Age = get_vals(files, z, axs, snapnum, i, on)
    print (met, Mratio, Type, Age)
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_user(met, Mratio, Type, z)

    x = 12 + np.log10(x)
    y = np.log10(y)
    xx = 12 + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    make_fig.fig_user(axs, z, x, y, xx, yy, yy_up, yy_low, den)

    xlim = [6.5,10.5]
    ylim = [-7,-0.1]
    xticks = [7, 8, 9, 10]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)

    return add


def plot_O_H_vs_DGR_median(files, z, axs, snapnum, i, on):

    add, met, Mratio, Type, Age = get_vals(files, z, axs, snapnum, i, on)

    xx, yy, yy_up, yy_low = create_out.out_median(met, Mratio, Type, z)

    xx = 12 + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    make_fig.fig_median(axs, z, xx, yy, yy_up, yy_low, i)

    xlim = [7,10.5]
    ylim = [-7,-0.1]
    xticks = [7, 8, 9, 10]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)

    return add


def plot_O_H_vs_DGR_age(files, z, axs, snapnum, i, on):

    add, met, Mratio, Type, Age = get_vals(files, z, axs, snapnum, i, on)

    x, y, xx, yy, yy_up, yy_low, den = create_out.out_age(met, Mratio, Type, Age, z)

    x = 12 + np.log10(x)
    y = np.log10(y)
    xx = 12 + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    p = make_fig.fig_age(axs, z, x, y, xx, yy, yy_up, yy_low, den)

    xlim = [7.5,10.5]
    ylim = [-7,-0.1]
    xticks = [8, 9, 10]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)

    return add, p, den
