import numpy as np
import matplotlib.pyplot as plt
filename = 'MagPhys3DHSTv01_best_agn0'
data = np.genfromtxt('./Obs_data/'+filename+'.csv', skip_header = 1, delimiter = ',', filling_values = np.nan, usecols = (1,2,3,4,5,6,7))

z, Mstar_l, Mstar, Mstar_u, Mdust_l, Mdust, Mdust_u = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6]
ok = np.logical_and(z >= 0., z <= 10)
z = z[ok]
Mstar = Mstar[ok]
Mdust = Mdust[ok]
Mstar_l = Mstar - Mstar_l[ok]
Mstar_u = Mstar_u[ok] - Mstar
Mdust_l = Mdust - Mdust_l[ok]
Mdust_u = Mdust_u[ok] - Mdust
z = np.round(z,0)
"""
bottom=0.09 
left = 0.08

fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize=(15, 13), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()
for z_ in range(0, 9):

    ok = np.logical_and(z == z_, z == z_)
    if np.sum(ok) >= 1:
        axs[z_].errorbar(Mstar[ok], Mdust[ok], xerr = [Mstar_l[ok], Mstar_u[ok]], yerr = [Mdust_l[ok], Mdust_u[ok]], fmt='.', color='indigo',label=r'$\mathrm{S.}$ $\mathrm{P.}$ $\mathrm{Driver}$ $\mathrm{2017}$', alpha = 0.05, elinewidth=0.1)
        axs[z_].grid()
        xlim = [7.5,11.9]
        ylim = [1.5,10.8]
        xticks = [8, 9, 10, 11]
        
        axs[z_].set_xlim(xlim)
        axs[z_].set_ylim(ylim)
        axs[z_].set_xticks(xticks)
fig.tight_layout()
fig.subplots_adjust(bottom=bottom, left = left, wspace=0, hspace=0)

plt.show()
"""
