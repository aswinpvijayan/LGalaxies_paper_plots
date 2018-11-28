#!/bin/bash 

# -- The shell used to interpret this script
#$ -S /bin/bash
# Do not lmit the stacksize (the maximum memory a job can use)
ulimit -s unlimited
# -- Execute this job from the current working directory.
#$ -cwd
#$ -pe openmpi 4
# Say which queue you want to submit to
#$ -q mps.q@@compute_amd
# Avoid node 106, 105, 104
#$ -l h=!(node103|node106|node105)
## Define a task farm of jobs
##$ -t 2-4
# Limit to 4 concurrent jobs
##$ -tc 4
# -- Job output to stderr will be merged into standard out. Remove this line if
# -- you want to have separate stderr and stdout log files
#$ -j y
#$ -o output/
# Give the job a name
#$ -N dust_prod_plot
# -- Send email when the job exits, is aborted or suspended
# #$ -m eas
# #$ -M YOUR_USERNAME@sussex.ac.uk


module load mps/software/

echo "Running job script"


#i=$(($SGE_TASK_ID))

#python phase_space_plots_user.py 4 5
#python DMF.py 5
python dust_prod_rates.py 1 0 5

echo "Finished job script"
