
import sys
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
from astropy.cosmology import Planck15
sys.path.append('/lustre/scratch/astro/ap629/my_lib/')
from essent import *

def outputs(x, y, sSFR, Type, z):
    
    """
    #Function remove_ is for removing any nan's, inifnities and other non-physical values, can specify if the zero values in the array should not be removed, by giving the array number as argument. It then selects the central halos which are above a particular sSFR. Returns the median (in 13 bins in the x axis in logarithmic space) with the values of the 16th percentile and the 84th percentile.
    """
    
    out = remove_(np.array([x, y, sSFR, Type]), np.array([3]))  
    del x, y, sSFR, Type
    out = out[(out[2] > sSFR_cut(z)) & (out[3] == 0)]
    out = np.array(out).T
    thisx = out[0]
    thisy = out[1]
    
    del out
    
    xx, yy, yy_up, yy_low = get_median(thisx, thisy, n = 20)
    
    return xx, yy, yy_up, yy_low

h = 0.673 #little h as used in the model

inp = int(input('Enter the appropriate number for the required phase space plot: \n0 - Stellar Mass versus Dust Mass \n1 - Metallicity versus Dust Mass\n2 - Dust fraction in cold gas versus Stellar Mass\n3 - Dust fraction in cold gas versus Metallicity\n'))

norm_z = np.arange(0, 4, 1)
high_z = np.arange(9, 14, 1)
all_z = np.arange(0, 14, 1)

if inp == 0:
    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(M_{Dust}/(M_{\odot}))$'
    savename = 'Dust_Stellar_MR'
elif inp == 1:
    xlab = r'$12+log_{10}(O/H)$'
    ylab = r'$log_{10}(M_{Dust}/(M_{\odot}))$'
    savename = 'Dust_metal_MR'
elif inp == 2:
    xlab = r'$log_{10}(M_{*}/(M_{\odot}))$'
    ylab = r'$log_{10}(M_{Dust}/(M_{cold\mathrm{ }gas}))$'
    savename = 'Dust_gas_ratio_MR.png'
elif inp == 3:
    xlab = r'$12+log_{10}(O/H)$'
    ylab = r'$log_{10}(M_{Dust}/(M_{cold\mathrm{ }gas}))$'
    savename = 'Dust_gas_ratio_MR'
else:
    print ('Not an applicable choice, retry....')
    sys.exit()  

z_inp = int(input('Select redshift (z) range: \n0 - [0, 8]\n1 - [9, 13]\n2 - [0, 13]\n'))

if z_inp == 0:
    z_ = norm_z
    name_z = 'z0_8'
elif z_inp == 1:
    z_ = high_z
    name_z = 'z9_13'
elif z_inp == 2:
    z_ = all_z
    name_z = 'z0_13'
else:
    print ('Not an applicable choice, retry....')
    sys.exit()  
""" 
turn_on = int(input('Input 1 in order to include MRII:\n'))
if turn_on == 1:
    name_x = '_MRII'
else:
    name_x = ''
"""
name_x = ''
Vars = np.array([['StellarMass', 'Dust_elements'], ['ColdGas_elements', 'Dust_elements'], ['ColdGas_elements', 'Dust_elements'], ['ColdGas_elements', 'Dust_elements']])
colours = ['blue', 'green', 'red', 'magenta', 'brown', 'orange', 'violet', 'cyan', 'indigo', 'black', 'yellow', 'blue', 'green', 'red']


fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(9, 8), sharex=True, sharey=True, facecolor='w', edgecolor='k')

for i, z in enumerate(z_):
    
    if turn_on == 1:
        file_req = '/lustre/scratch/astro/ap629/Dust_output_17jan/MR/hdf5/lgal_z{}_N*.hdf5'.format(z)
    else:
        file_req = '/research/astro/virgo/Aswin/data/MR/hdf5/lgal_z{}_N*.hdf5'.format(z)
    
    Mstar = get_data1('StellarMass', file_req)*(1e10)/h  #Now in units of Msun 
    SFR = get_data1('Sfr', file_req)    #Star formation rate, in units of Msun/yr
    sSFR = SFR/Mstar    #Specific star formation rate, in units of /yr
    del SFR
    Type = get_data1('Type', file_req)
    
    if inp == 0:
        
        Mdust = get_data1(Vars[inp][1], file_req)  #Array of mass of 11 elements in dust: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
        Mdust = np.sum(Mdust, axis = 1)  #Total dust mass
        
        xx, yy, yy_up, yy_low = outputs(Mstar, Mdust, sSFR, Type, z)
        
        xx = np.log10(xx)
        yy = np.log10(yy)
        yy_up = np.log10(yy_up)
        yy_low = np.log10(yy_low)
        
    if inp == 1:
        
        Mcg = get_data1(Vars[inp][0], file_req)  #Array of mass of 11 elements present in the cold gas: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
        Mdust = get_data1(Vars[inp][1], file_req) #Array of mass of 11 elements in dust: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
        Mdust = np.sum(Mdust, axis = 1)  #Total dust mass
        met = Mcg[:,4]/(16.0*Mcg[:,0])  #Metallicity = Number of Oxygen atoms/Number of Hydrogen atoms
        xx, yy, yy_up, yy_low = outputs(met, Mdust, sSFR, Type, z)
        
        xx = 12.0 + np.log10(xx)
        yy = np.log10(yy)
        yy_up = np.log10(yy_up)
        yy_low = np.log10(yy_low)
        
    if inp == 2:
        
        Mcg = get_data1(Vars[inp][0], file_req) #Array of mass of 11 elements present in the cold gas: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
        Mcg = np.sum(Mcg, axis = 1) #Total cold gas mass
        Mdust = get_data1(Vars[inp][1], file_req) #Array of mass of 11 elements in dust: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
        Mdust = np.sum(Mdust, axis = 1) #Total dust mass
        Mratio = Mdust/Mcg   #Dust-to-total cold gas mass ratio. Dimensionless.
        xx, yy, yy_up, yy_low = outputs(Mstar, Mratio, sSFR, Type, z)
        
        xx = np.log10(xx)
        yy = np.log10(yy)
        yy_up = np.log10(yy_up)
        yy_low = np.log10(yy_low)
        
    if inp == 3:
        
        Mcg = get_data1(Vars[inp][0], file_req) #Array of mass of 11 elements present in the cold gas: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
        met = Mcg[:,4]/(16.0*Mcg[:,0]) #Metallicity = Number of Oxygen atoms/Number of Hydrogen atoms
        Mcg = np.sum(Mcg, axis = 1) #Total cold gas mass
        Mdust = get_data1(Vars[inp][1], file_req) #Array of mass of 11 elements in dust: H, He, C, N, O, Ne, Mg, Si, S, Ca and Fe. In units of Msun
        Mdust = np.sum(Mdust, axis = 1)  #Total dust mass
        Mratio = Mdust/Mcg  #Dust-to-total cold gas mass ratio. Dimensionless.
        xx, yy, yy_up, yy_low = outputs(met, Mratio, sSFR, Type, z)
        
        xx = 12.0 + np.log10(xx)
        yy = np.log10(yy)
        yy_up = np.log10(yy_up)
        yy_low = np.log10(yy_low)
        
    axs.plot(xx, yy, lw = 2, label = r'$z = {}$'.format(z), color = colours[i])
    axs.plot(xx, yy_up, lw = 2, ls = 'dashed', color = colours[i])
    axs.plot(xx, yy_low, lw = 2, ls = 'dashed', color = colours[i])
    
    for label in (axs.get_xticklabels() + axs.get_yticklabels()):
        label.set_fontsize(18)
    
    axs.grid()
    lgd = axs.legend(frameon=False, fontsize = 18, markerscale=2, loc = 4, numpoints=1, handletextpad=0.005)

    
axs.set_xlabel(xlab, fontsize=22)
axs.set_ylabel(ylab, fontsize=22)
fig.subplots_adjust(bottom=0.1, left = 0.1, wspace=0, hspace=0)
fig.tight_layout()
plt.savefig(name_z+savename+name_x+'median.png')
plt.show()
