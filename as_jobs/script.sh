#!/bin/bash 
######################################################################
# Options for the batch system
# These options are not executed by the script, but are instead read by the
# batch system before submitting the job. Each option is preceeded by '#$' to
# signify that it is for grid engine.
#
# All of these options are the same as flags you can pass to qsub on the
# command line and can be **overriden** on the command line. see man qsub for
# all the details
######################################################################
# -- The shell used to interpret this script
#$ -S /bin/bash
# Do not lmit the stacksize (the maximum memory a job can use)
ulimit -s unlimited
# -- Execute this job from the current working directory.
#$ -cwd
# Say which queue you want to submit to
#$ -q mps.q@@compute_amd
# Define a task farm of jobs
#$ -t 1-4
# Limit to 4 concurrent jobs
#$ -tc 4
# -- Job output to stderr will be merged into standard out. Remove this line if
# -- you want to have separate stderr and stdout log files
#$ -j y
#$ -o output/
# Give the job a name
#$ -N dust_plots
# -- Send email when the job exits, is aborted or suspended
# #$ -m eas
# #$ -M YOUR_USERNAME@sussex.ac.uk

######################################################################
# Job Script
# Here we are writing in bash (as we set bash as our shell above). In here you
# should set up the environment for your program, copy around any data that
# needs to be copied, and then execute the program
######################################################################
module load mps/software/
# Here we execute usual shell commands like any other shell script. The
# output here will be sent to the job's standard out
echo "Running job script"

# Finally we run our executable. Here we are passing the command line argument
# above to the script

i=$(($SGE_TASK_ID -1))

python dust_prod_rates.py $i $i 0
#python dust_prod_rates.py $i $i 1

python phase_space_plots.py $i 0
#python phase_space_plots.py $i 1

#python phase_space_plots_median_zframes.py $i 

python phase_space_plots_median.py $i 1 0

#python phase_space_plots_median.py $i 0 1

#python phase_space_plots_median.py $i 1 1

#if [ $i -lt 3 ]
#then
    python dust_prod_rates.py $i 0 0
#    python dust_prod_rates.py $i $i 1
    python props_redshift.py $i 0
#    python props_redshift.py $i 1
#fi


#if [ $i -lt 2 ]
#then
    python DMF.py $i
#fi
#
echo "Finished job script"
