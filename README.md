# LGalaxies_Dust paper plots

This repository contains the scripts for analysing the output of the L-Galaxies SAMs implementation of dust production in its galaxies. The code producing the data can be found at https://github.com/aswinpvijayan/L-Galaxies_Dust which uses the dark matter merger trees of the Millennium (MR) and the Millennium II (MRII) runs.

* All the masses used in the script are in Msolar. Hubble parameter, h = 0.673 is used since that is the value used in the SAMs. Most of the plots use both the MR and MRII data, hence the data collection routine uses 4 cpus to read data using the joblib python module.

* Figure 1 is produced using frac_plots.py
* Figures 2, 4, 5, 6, 7 and A1 produced using phase_space_plots_user.py
* Figure 3 produced using phase_space_plots_age.py
* Figures 8, 9 and 13c produced using dust_prod_rates.py
* Figure 10 produced using rates_compare.py
* Figures 11, 13a and 13b produced using phase_space_plots_median.py
* Figure 12 produced using DMF.py
* Figure A2 produced using diff_taccs.py
