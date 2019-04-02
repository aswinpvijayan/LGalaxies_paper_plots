# LGalaxies_Dust paper plots

This repository contains the scripts for analysing the output of the L-Galaxies SAMs implementation of dust production in its galaxies. The code producing the data can be found at https://github.com/aswinpvijayan/L-Galaxies_Dust which uses the dark matter merger trees of the Millennium (MR) and the Millennium II (MRII) runs. When running the scripts change the ```fileMR``` and ```fileMRII``` to the respective locations of the output. Folder 'func_def' contains all the required python files that are being called from the main directory.

* All the masses used in the script are in Msolar. Hubble parameter, h = 0.673 is used since that is the value used in the SAMs. Most of the plots use both the MR and MRII data, hence the data collection routine uses 4 cpus to read data using the joblib python module.

* Figure 1 is produced using frac_plots.py
* Figure 2 is produced using depletions.py
* Figures 3, 6, 10, 11, 12 and A1 produced using phase_space_plots_user.py
* Figure 4 produced using DTM_1d_hist.py
* Figure 5 produced using phase_space_plots_age.py
* Figure 7 produced using scaling_DTM.py
* Figures 8 and 9 produced using dust_prod_rates.py
* Figure 13 produced using DMF.py
* Figure A2 produced using diff_taccs.py
* Figure A3 is produced using fit_plot.py
