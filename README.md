# LGalaxies_Dust paper plots

This repository contains the scripts for analysing the output of the L-Galaxies SAMs implementation of dust production in its galaxies. The code producing the data can be found at https://github.com/aswinpvijayan/L-Galaxies_Dust which uses the dark matter merger trees of the Millennium (MR) and the Millennium II (MRII) runs. When running the scripts change the ```fileMR``` and ```fileMRII``` to the respective locations of the output.

* All the masses used in the script are in Msolar. Hubble parameter, h = 0.673 is used since that is the value used in the SAMs. Most of the plots use both the MR and MRII data, hence the data collection routine uses 4 cpus to read data using the joblib python module.

* Figure 1 is produced using frac_plots.py
* Figures 2, 4, 8, 9, 10 and A1 produced using phase_space_plots_user.py
* Figure 3 produced using phase_space_plots_age.py
* Figures 6, 7 and 12c produced using dust_prod_rates.py
* Figure 11 produced using DMF.py
* Figures 12a and 12b produced using phase_space_plots_median.py
* Figure A2 produced using diff_taccs.py
