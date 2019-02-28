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
        ok = np.where(Mstar >= 10**9.0)[0]
    elif on and i == 1:
        ok = np.logical_and(Mstar > 10**7.0, Mstar < 10**9.0)
    else:
        ok = np.array([True]*len(Mstar))

    Mstar = Mstar[ok]
    Type = get_.get_var(files[i], 'Type', snap)[ok]
    Age = get_.get_var(files[i], 'MassWeightAge', snap)[ok]
    Mcg = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok] + get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]
    Mmet = np.nansum(Mcg[:,2:], axis = 1)
    Mcg = np.nansum(Mcg, axis = 1)
    Z = Mmet/Mcg
    tacc = get_.get_var(files[i], 't_acc', snap)[ok]
    
    ok = np.where(Mcg > 1e6)  #Only selecting galaxies with cold gas mass above this value, so the galaxies we are looking at are realistic for the dust content usually seen
    
    return Mstar[ok], tacc[ok], Type[ok], Age[ok], Z[ok]


def plot_tacc_Mstar_user(files, z, axs, snapnum, i, on):

    add = sims[i]

    Mstar, tacc, Type, Age, Z = get_vals(files, z, axs, snapnum, i, on)
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_age(Mstar, tacc, Type, Z, z)

    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)
    den = np.log10(den)
    
    #p = make_fig.fig_age(axs, z, x, y, xx, yy, yy_up, yy_low, den, bins = None) #Colour by metallicity
    p = make_fig.fig_user(axs, z, x, y, xx, yy, yy_up, yy_low, den)
    
    xlim = [7.5,11.7]
    ylim = [5,11.5]
    xticks = [8, 9, 10, 11]
    yticks = [5, 6, 7, 8, 9, 10, 11]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    axs.set_yticks(yticks)

    return add, p, den


def plot_tacc_Mstar_median(files, z, axs, snapnum, i, on):

    add = sims[i]

    Mstar, tacc, Type, Age, Z = get_vals(files, z, axs, snapnum, i, on)
    xx, yy, yy_up, yy_low = create_out.out_median(Mstar, tacc, Type, z)

    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    make_fig.fig_median(axs, z, xx, yy, yy_up, yy_low, i)

    xlim = [7.5,11.7]
    ylim = [5,10.5]
    xticks = [8, 9, 10, 11]
    yticks = [5, 6, 7, 8, 9, 10]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    axs.set_yticks(yticks)

    return add


def plot_tacc_Mstar_age(files, z, axs, snapnum, i, on):

    add = sims[i]

    Mstar, tacc, Type, Age, Z = get_vals(files, z, axs, snapnum, i, on)
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_age(Mstar, tacc, Type, Age, z)

    x = np.log10(x)
    y = np.log10(y)
    xx = np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    p = make_fig.fig_age(axs, z, x, y, xx, yy, yy_up, yy_low, den)

    xlim = [7.5,11.7]
    ylim = [5,10.5]
    xticks = [8, 9, 10, 11]
    yticks = [5, 6, 7, 8, 9, 10]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)

    return add, p, den
