import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
sys.path.append('../')
import matplotlib as mpl
import matplotlib.pyplot as plt


def fig_user(axs, z, x, y, xx, yy, yy_up, yy_low, den):
    
    gridsize = np.array([(int(max(x)-min(x))/0.1), int((max(y)-min(y))/0.1)]).astype(int)
    hb = axs.hexbin(x, y, gridsize=gridsize, bins = 'log', cmap = plt.cm.get_cmap('gist_yarg'), mincnt = 1)
    #mini = hb.norm.vmin
    #maxi = hb.norm.vmax
    #normalize = mpl.colors.Normalize(vmin=mini, vmax=maxi)
    #axs.hexbin(x, y, gridsize=gridsize, bins = 'log', cmap = plt.cm.get_cmap('gist_yarg'), mincnt = 1, norm = normalize)
    
    axs.plot(xx, yy, lw = 1, color = 'orange')
    axs.plot(xx, yy_up, lw = 1, ls = 'dashed', color = 'orange')
    axs.plot(xx, yy_low, lw = 1, ls = 'dashed', color = 'orange')
    
    del x, y, xx, yy, yy_up, yy_low
    
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(16)
    
    axs.grid(True)
    lgd = axs.legend(fontsize = 15, loc = 0, markerscale=2, numpoints=1, handletextpad=0.005)
    if np.isscalar(lgd):
        lgd.set_zorder(100)
        
        
def fig_median(axs, z, xx, yy, yy_up, yy_low, i):
    
    colours = ['blue', 'green', 'red', 'magenta', 'brown', 'orange', 'violet', 'cyan', 'indigo', 'black', 'violet', 'blue', 'green', 'red', 'magenta', 'brown']
    
    if i == 1:
        label = r'$z={}$'.format(z)
    else:
        label='_nolegend_'
        
    axs.errorbar(xx, yy, yerr = [yy - yy_low, yy_up - yy], lw = 2, color = colours[z], label = label)
    
    try:
        lgd = axs.legend(fontsize = 15, markerscale=2, loc = 4, numpoints=1, handletextpad=0.005, frameon = False)
        if np.isscalar(lgd):
            lgd.set_zorder(100)
    except:
        pass
    
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(16)
    
    axs.grid(True)
    
    
def fig_age(axs, z, x, y, xx, yy, yy_up, yy_low, den, bins = 'log'):
    
    gridsize = np.array([(int(max(x)-min(x))/0.02), int((max(y)-min(y))/0.02)]).astype(int)
    p = axs.hexbin(x, y, gridsize=gridsize, bins = bins, alpha = 0.7, C = den, cmap = plt.cm.get_cmap('jet'), mincnt = 1)
    
    axs.plot(xx, yy, lw = 1, color = 'brown')
    axs.plot(xx, yy_up, lw = 1, ls = 'dashed', color = 'brown')
    axs.plot(xx, yy_low, lw = 1, ls = 'dashed', color = 'brown')
    
    del x, y, xx, yy, yy_up, yy_low
    
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(16)
    
    axs.grid(True)
    lgd = axs.legend(fontsize = 15, markerscale=2, loc = 4, numpoints=1, handletextpad=0.005)
    if np.isscalar(lgd):
        lgd.set_zorder(100)

    return p

