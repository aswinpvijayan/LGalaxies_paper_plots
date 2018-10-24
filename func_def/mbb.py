import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants
from scipy import integrate

_h = constants.h.si
_c = constants.c.si
_k_B = constants.k_B.si

def _blackbody_hz(nu, temperature):
    """
    Compute the Planck function given nu in Hz and temperature in K with output
    in cgs
    """
    I = (2*_h*(nu*u.Hz)**3 / _c**2) * (np.exp(_h*nu*u.Hz/(_k_B*temperature*u.K)) - 1)**-1

    return I
    
def Fac(nu):

    return (0.0383*(nu/(2.216*1e12))**2) * (u.m**2 / u.kg)
    

def _MBB(nu, temperature):

    return Fac(nu)*_blackbody_hz(nu, temperature)
    

def prefactor(temperature):
    
    
    
    nu1 = 6.279*1e14
    nu2 = 1.884*1e12

    integral = lambda x, a : _MBB(x, a).si.value

    y, err = integrate.quad(integral, nu2, nu1, args = (temperature, ))
    y = y/(3.848*1e26)
    #y = y*(u.kg)*(u.s**(-3)) 

    y = y * 1.988435*1e30
    return y
#Factor = 887.058*4pi
