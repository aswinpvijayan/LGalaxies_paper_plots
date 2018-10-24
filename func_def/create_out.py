import sys
import os
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
sys.path.append('../')
import get_

"""
    Function remove_ is for removing any nan's, inifnities and other non-physical values, can specify if the zero values in the array should not be removed, by giving the array number as argument. It then selects the central halos. 
"""

def out_user(x, y, Type, z):
    
    out = get_.remove_(np.array([x, y, Type]), np.array([2]))  
    out = out[(out[2] == 0)]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    
    del out, x, y, Type
    
    xx, yy, yy_up, yy_low = get_.get_median(thisx, thisy, n = 10)
    
    den = get_.get_density(thisx, thisy)  #Returns the number density in a pixel and the median with the values of the 16th percentile and the 84th percentile. 
                                          
    return thisx, thisy, xx, yy, yy_up, yy_low, den
    
    
def out_median(x, y, Type, z):
    
    out = get_.remove_(np.array([x, y, Type]), np.array([2]))  
    out = out[(out[2] == 0)]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    
    del out, x, y, Type
    
    xx, yy, yy_up, yy_low = get_.get_median(thisx, thisy, n = 15)
    
    return xx, yy, yy_up, yy_low
    

def out_age(x, y, Type, Age, z):

    out = get_.remove_(np.array([x, y, Age, Type]), np.array([3]))  
    out = out[(out[2] > 0) & (out[3] == 0)]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    den = out[2]
    
    del out, x, y, Age, Type
    
    xx, yy, yy_up, yy_low = get_.get_median(thisx, thisy, n = 10)
    
    return thisx, thisy, xx, yy, yy_up, yy_low, den
