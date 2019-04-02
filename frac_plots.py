#!/usr/bin/env python

"""

    Figure 1 in paper

"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("paper")
sns.set_style(style='white') 
# t_e
t_e_array = np.array([0.01, 1.0, 100])

# mu
mu_array = np.array([0.05,0.9])

# Initial dust fractions
f0_array = np.array([0.05,0.05])

# Time range for integration
t=np.logspace(-3,3,4000)

# Plotting
#sns.set_context('paper')


# Define the derivative function
def dfdt(f,t,t_e,mu,fcmax):
    fc,fi=f
    return [fcmax-fc+(fi-fc)/t_e, mu*(fc-fi)/((1-mu)*t_e)]


fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(13, 8), sharex=True, sharey=True, facecolor='w', edgecolor='k')

f0=np.array([f0_array[0],f0_array[1]])

for i in range(0, len(mu_array)):

    mu = mu_array[i]
    for j in range(0, len(t_e_array)):

        t_e = t_e_array[j]

        f=odeint(dfdt,f0,t,args=(t_e,mu, 0.7))
        fc=f[:,0]
        fi=f[:,1]
        
        if i == 0 and j == 0:
            axs[i, j].semilogx(t,fc,color='blue',label=r'$\mathrm{f_{C}}$')
            axs[i, j].semilogx(t,fi,color='green',label=r'$\mathrm{g_C}$')
            axs[i, j].plot([t_e,t_e],[0.,1.],':',color='red')
            axs[i, j].plot([t_e*(1.-mu)/mu,t_e*(1.-mu)/mu],[0.,1.],'-.',color='orange')
            axs[i, j].axhline(y = 0.7, ls = '-.',label = r'$\mathrm{f_{C,max}}$',color='blue')
            axs[i, j].legend(frameon=True, fontsize = 19, loc = 0, framealpha = 0.7)
        if i == 0 and j == 1:
            axs[i, j].semilogx(t,fc,color='blue')
            axs[i, j].semilogx(t,fi,color='green')
            axs[i, j].plot([t_e,t_e],[0.,1.],':',color='red',label=r'$\frac{\tau_{\mathrm{exch}}}{\mathrm{t_{acc}}}$')
            axs[i, j].plot([t_e*(1.-mu)/mu,t_e*(1.-mu)/mu],[0.,1.],'-.',color='orange',label=r'$\frac{\tau_{\mathrm{exch}}^\prime}{\mathrm{t_{acc}}}$')
            axs[i, j].axhline(y = 0.7, ls = '-.',color='blue')
            axs[i, j].legend(frameon=True, fontsize = 25, loc = 0, framealpha = 0.7)
        else:
            axs[i, j].semilogx(t,fc,color='blue')
            axs[i, j].semilogx(t,fi,color='green')
            axs[i, j].plot([t_e,t_e],[0.,1.],':',color='red')
            axs[i, j].plot([t_e*(1.-mu)/mu,t_e*(1.-mu)/mu],[0.,1.],'-.',color='orange')
            axs[i, j].axhline(y = 0.7, ls = '-.',color='blue')
        #axs.set_xlabel(r'$t/t_{acc}$', fontsize = 20)
        #axs.set_ylabel(r'$f_C,\ f_D$', fontsize = 20)
        for label in (axs[i, j].get_xticklabels() + axs[i, j].get_yticklabels()):
            label.set_fontsize(18)

        axs[i, j].set_xlim((7e-3,7e2))
        axs[i, j].set_ylim((0,1))
        axs[i, j].set_yticks(np.arange(0,1.2,0.2))
        #axs[i, j].grid(True)
        
            


fig.subplots_adjust(bottom=0.1, left = 0.1, wspace=0, hspace=0.13)

fig.text(0.12, 0.92, r'$\mu=0.05,\ \tau_{\mathrm{exch}}/\tau_{\mathrm{acc}}=0.01$', va='center', fontsize=19)
fig.text(0.40, 0.92, r'$\mu=0.05,\ \tau_{\mathrm{exch}}/\tau_{\mathrm{acc}}=1$', va='center', fontsize=19)
fig.text(0.67, 0.92, r'$\mu=0.05,\ \tau_{\mathrm{exch}}/\tau_{\mathrm{acc}}=100$', va='center', fontsize=19)

fig.text(0.12, 0.495, r'$\mu=0.9,\ \tau_{\mathrm{exch}}/\tau_{\mathrm{acc}}=0.01$', va='center', fontsize=19)
fig.text(0.40, 0.495, r'$\mu=0.9,\ \tau_{\mathrm{exch}}/\tau_{\mathrm{acc}}=1$', va='center', fontsize=19)
fig.text(0.67, 0.495, r'$\mu=0.9,\ \tau_{\mathrm{exch}}/\tau_{\mathrm{acc}}=100$', va='center', fontsize=19)

fig.text(0.03, 0.5, r'$\mathrm{f_{C}},\ \mathrm{g_{C}}$', va='center', rotation='vertical', fontsize=22)
fig.text(0.5, 0.035, r'$\mathrm{t}/\tau_{\mathrm{acc}}$', va='center', fontsize=22)

plt.savefig('fc_gc.pdf', bbox_inches='tight')
plt.show()
