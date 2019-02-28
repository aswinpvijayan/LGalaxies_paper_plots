import glob
import os.path
import re
import numpy as np
import pandas as pd
import h5py
from joblib import Parallel, delayed
import timeit


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):

    return [atoi(c) for c in re.split('(\d+)', text)]

def get_files(file_req):
    
    if np.isscalar(file_req):
        files = glob.glob(file_req)
        files.sort(key=natural_keys)
    else:
        files = []
        for i in file_req:
            tmp = glob.glob(i)
            tmp.sort(key=natural_keys)
            files.append(tmp)
        files = np.concatenate(files)
    return files
    
def get_(filename, pack):
    
    var, snap = pack
    #print ('Reading in {} from {}...'.format(var, filename))
    f_0 = h5py.File(filename, 'r')
    g = f_0[snap]
    out = np.array(g[var])
    #print (len(out))
    
        
    return out


def get_var(files, var, snap):

    start = timeit.default_timer()
    print ('Reading in {} from snap {}...'.format(var, snap))    
    filess = get_files(files)
    
    if snap in ['17', '16', '15', '14', '13', '12', '11', '10', '9', '8', '7']:
        
        kill = np.array([])
        out = Parallel(n_jobs = 4)(delayed(get_)(i, ['Type', snap]) for i in filess)
        
        for i in range(0, len(out)):
            if len(out[i]) == 0:
                
                kill = np.append(kill, i)
        
        if len(kill) != 0:
            
            filess = np.delete(filess, kill)
        
    out = Parallel(n_jobs = 16)(delayed(get_)(i, [var, snap]) for i in filess)  
    out = np.concatenate(out, axis = 0)
    
    stop = timeit.default_timer()

    print ('Time taken: {}'.format(stop - start))
    
    return out
    
    
def remove_(x, turn_off = [10**20]):
    
    """
    A function to remove the non-essential elements in related datasets, with option for 
    turning off if a certain dataset's zero values should be removed or not.
    """
    
    df = pd.DataFrame(data = x.T)
    
    for i in range(0, len(df.columns)):
        ok = np.logical_not(np.isinf(df[i]))
        df = df.ix[ok]
        df = pd.DataFrame(data = np.array(df), index = range(0,len(df[0])))
        ok = np.logical_not(np.isnan(df[i]))
        df = df.ix[ok]
        df = pd.DataFrame(data = np.array(df), index = range(0,len(df[0])))
    
        if i not in turn_off:
            ok = np.where(df[i] > 0)[0]
            df = df.ix[ok]
            df = pd.DataFrame(data = np.array(df), index = range(0,len(df[0])))
    
    return df
    
def get_median(x, y, n = 15):
    
    #bins = 10**(np.arange(min(np.log10(x))-0.2, max(np.log10(x))+0.2, 0.3))
    bins = 10**(np.linspace(min(np.log10(x)), max(np.log10(x)), num = n))
    xx = yy = yy_up = yy_low = np.array([])

    for i in range(0, len(bins)-1):
        ok = np.logical_and(x>=bins[i], x<bins[i+1])
        if np.sum(ok)>5:
            xx = np.append(xx, np.nanmedian(x[ok]))
            yy = np.append(yy, np.nanmedian(y[ok]))
            
            yy_low = np.append(yy_low, np.percentile(y[ok], 16))
            yy_up = np.append(yy_up, np.percentile(y[ok], 84))

    return xx, yy, yy_up, yy_low

def get_density(x, y):
    
    x_edges = 10**(np.arange(min(np.log10(x))-0.1, max(np.log10(x))+0.1, 0.05))
    y_edges = 10**(np.arange(min(np.log10(y))-0.1, max(np.log10(y))+0.1, 0.05))
    
    hist, xedges, yedges = np.histogram2d(x, y, bins=(x_edges, y_edges))
    xidx = np.digitize(x, x_edges)-1
    yidx = np.digitize(y, y_edges)-1
    norm = (np.sum(hist**2))**(0.5)
    num_pts = (hist[xidx, yidx])/norm
    
    return num_pts

from astropy.cosmology import Planck15
from astropy import units as u

def sSFR_cut(z):
        
    tmp = 1/Planck15.H(z).decompose()
    return (1.0/(3*tmp.to(u.yr))).value
    
