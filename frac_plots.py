import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns

# t_e
t_e_array = np.array([0.01, 1.0, 100])

# mu
mu_array = np.array([0.05,0.9])

# Initial dust fractions
f0_array = np.array([0.05,0.05])

# Time range for integration
t=np.logspace(-3,3,4000)

# Plotting
sns.set_context('paper')


# Define the derivative function
def dfdt(f,t,t_e,mu,fcmax):
    fc,fi=f
    return [fcmax-fc+(fi-fc)/t_e, mu*(fc-fi)/((1-mu)*t_e)]
    
    
fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(12, 7), sharex=True, sharey=True, facecolor='w', edgecolor='k')

f0=np.array([f0_array[0],f0_array[1]])

for i in range(0, len(mu_array)):
    
    mu = mu_array[i]
    for j in range(0, len(t_e_array)):

        t_e = t_e_array[j]

        f=odeint(dfdt,f0,t,args=(t_e,mu, 0.3))
        fc=f[:,0]
        fi=f[:,1]

        axs[i, j].semilogx(t,fc,label=r'$f_C$')
        axs[i, j].semilogx(t,fi,label=r'$f_D$')
        axs[i, j].plot([t_e,t_e],[0.,1.],':',label=r'$\frac{\tau_{exch}}{t_{acc}}$')
        axs[i, j].plot([t_e*(1.-mu)/mu,t_e*(1.-mu)/mu],[0.,1.],'-.',label=r'$\frac{\tau_{exch}^\prime}{t_{acc}}$')
        axs[i, j].axhline(y = 0.2, ls = '-.', label = r'$f_{C,max}$')

        #axs.set_xlabel(r'$t/t_{acc}$', fontsize = 20)
        #axs.set_ylabel(r'$f_C,\ f_D$', fontsize = 20)
        for label in (axs[i, j].get_xticklabels() + axs[i, j].get_yticklabels()):
            label.set_fontsize(18)
        
        axs[i, j].set_xlim((7e-3,7e2))
        axs[i, j].set_ylim((0,1.1))
        axs[i, j].set_yticks(np.arange(0,1.2,0.2))
        axs[i, j].grid(True)
        if i == 0 and j == 0:
            axs[i, j].legend(frameon=False, fontsize = 16, loc = 0)

   
fig.subplots_adjust(bottom=0.1, left = 0.1, wspace=0, hspace=0)
fig.text(0.01, 0.5, r'$f_C,\ f_D$', va='center', rotation='vertical', fontsize=22)
fig.text(0.5, 0.035, r'$t/t_{acc}$', va='center', fontsize=22)

plt.savefig('fc_fd.eps')        
plt.show()
