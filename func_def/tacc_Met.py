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
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        ok = np.where(Mstar >= 10**8.85)[0]
        Mstar = Mstar[ok]
        Type = get_.get_var(files[i], 'Type', snap)[ok]
        Age = get_.get_var(files[i], 'MassWeightAge', snap)[ok]
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok]
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]
        Mcg = Mcg1 + Mcg2

        met = (Mcg[:,4])/(15.9994*(Mcg[:,0]))  #O/H ratio

        tacc = get_.get_var(files[i], 't_acc', snap)[ok]

        i = 1
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        tmp = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673
        ok = np.logical_and(tmp > 10**7.5, tmp < 10**8.9)
        Mstar = np.append(Mstar, tmp[ok])
        Type = np.append(Type, get_.get_var(files[i], 'Type', snap)[ok])
        Age = np.append(Age, get_.get_var(files[i], 'MassWeightAge', snap)[ok])
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)[ok]
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)[ok]

        met = np.append(met, (Mcg1[:,4]+Mcg2[:,4])/(15.9994*(Mcg1[:,0]+Mcg2[:,0])))

        tacc = np.append(tacc, get_.get_var(files[i], 't_acc', snap)[ok])

        add = 'MR_MRII'
    else:
        add = sims[i]
        snap = snapnum[i][np.where(redshift == str(z))[0][0]]
        Mstar = (get_.get_var(files[i], 'StellarMass', snap)*1e10)/0.673

        Type = get_.get_var(files[i], 'Type', snap)
        Mcg1 = get_.get_var(files[i], 'ColdGasDiff_elements', snap)
        Mcg2 = get_.get_var(files[i], 'ColdGasClouds_elements', snap)
        Mcg = Mcg1 + Mcg2

        met = (Mcg[:,4])/(15.9994*(Mcg[:,0]))  #O/H ratio

        tacc = get_.get_var(files[i], 't_acc', snap)

    return add, met, tacc, Type, Age


def plot_tacc_Met_user(files, z, axs, snapnum, i, on):

    add, met, tacc, Type, Age = get_vals(files, z, axs, snapnum, i, on)
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_user(met, tacc, Type, z)

    x = 12. + np.log10(x)
    y = np.log10(y)
    xx = 12. + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    make_fig.fig_user(axs, z, x, y, xx, yy, yy_up, yy_low, den)

    xlim = [7.5,10.5]
    ylim = [5,10.5]
    xticks = [8, 9, 10]
    yticks = [5, 6, 7, 8, 9, 10]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    axs.set_yticks(yticks)

    return add


def plot_tacc_Met_median(files, z, axs, snapnum, i, on):

    add, met, tacc, Type, Age = get_vals(files, z, axs, snapnum, i, on)
    xx, yy, yy_up, yy_low = create_out.out_median(met, tacc, Type, z)

    xx = 12. + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    make_fig.fig_median(axs, z, xx, yy, yy_up, yy_low, i)

    xlim = [7.5,10.5]
    ylim = [5,10.5]
    xticks = [8, 9, 10]
    yticks = [5, 6, 7, 8, 9, 10]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)
    axs.set_yticks(yticks)

    return add


def plot_tacc_Met_age(files, z, axs, snapnum, i, on):

    add, met, tacc, Type, Age = get_vals(files, z, axs, snapnum, i, on)
    x, y, xx, yy, yy_up, yy_low, den = create_out.out_age(met, tacc, Type, Age, z)

    x = 12. + np.log10(x)
    y = np.log10(y)
    xx = 12. + np.log10(xx)
    yy = np.log10(yy)
    yy_up = np.log10(yy_up)
    yy_low = np.log10(yy_low)

    p = make_fig.fig_age(axs, z, x, y, xx, yy, yy_up, yy_low, den)

    xlim = [7.5,10.5]
    ylim = [5,10.5]
    xticks = [8, 9, 10]
    yticks = [5, 6, 7, 8, 9, 10]

    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xticks(xticks)

    return add, p, den
