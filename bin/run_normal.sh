#!/bin/bash
# Author: F.N. Krohg
# Last-updated: 28.05.19
# Slurm job script for launching chiral superconductivity calculation on the Fram
# supercomputer. Based on example MPI Job script

#################################
# Setting Slurm variables       #
#################################

## Required settings are --account, --time and --nodes

## Project:
#SBATCH --account=nn2819k

## Job name:
#SBATCH --job-name=LowField
## Allocating amount of resources:
#SBATCH --nodes=4
## Number of tasks (aka processes) to start on each node: Pure mpi, one task per core
#SBATCH --ntasks-per-node=32 --cpus-per-task=1
## No memory pr task since this option is turned off on Fram in QOS normal.
## Run for x minutes, syntax is d-hh:mm:ss
#SBATCH --time=7-00:00:00 

# turn on all mail notification, and also provide mail address:
##SBATCH --mail-type=ALL
##SBATCH -M fredrik.n.krohg@ntnu.no

# you may not place bash commands before the last SBATCH directive
######################################################
## Setting variables and prepare runtime environment:
##----------------------------------------------------
## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

# Needed for jobs in normal or optimist partition:
#export SLURM_MEM_PER_CPU=1920

#1-01:00:00 
error_exit()
{
	echo "$1" 1>&2
	exit 1
}

SOURCE_PATH="/cluster/projects/nn2819k/finite-temp-vortex-lattice/Source/"
SCRIPTS_PATH="/cluster/projects/nn2819k/finite-temp-vortex-lattice/Scripts/"
JULIA_PATH="/cluster/home/fredrkro/Programs/julia-1.0.4/bin/julia"

# Each node will run 3 different temperature intervals, where each temp. int. contains 3 different temperatures.
T_start1=(1.999 1.972 1.945 1.919) #(1.975) 
T_start2=(1.99 1.963 1.936 1.91)
T_start3=(1.981 1.954 1.927 1.901)
# Each interval will be initialized from a separate state with temperatures set below
T_init1=(2.0 2.0 2.0 2.0)
T_init2=(2.0 2.0 2.0 2.0)
T_init3=(2.0 2.0 2.0 2.0)


declare -a outputs
for (( i=0; i<${#T_start1[@]}; i++))
do
    outputs[i]="T_start1_${T_start1[i]}_field.out"
done

# Needs to be edited when changing scripts
JULIA_SCRIPT="driver_Ts.jl"
AUX_SCRIPT="script_functions.jl"
CPUS="50" # Max is --ntasks-per-node * --nodes * --cpus-per-task (is it though?)


# Making sure we are in the right directory
[ -d $SCRIPTS_PATH ] || error_exit "ERROR: Could not find Source directory in $SCRIPTS_PATH"
cd $SCRIPTS_PATH
[ -f ./$JULIA_SCRIPT ] || error_exit "ERROR: Could not find julia script $JULIA_SCRIPT in `pwd`"
cp $JULIA_SCRIPT $SCRATCH/
[ -f ./$AUX_SCRIPT ] || error_exit "ERROR: Could not find auxhillary script code $AUX_SCRIPT in `pwd`"
cp $AUX_SCRIPT $SCRATCH/

#######################################################
## Prepare jobs, moving input files and making sure 
# output is copied back and taken care of
##-----------------------------------------------------

# Prepare input files
#cp inputfiles $SCRATCH
#cd $SCRATCH

# Make sure output is copied back after job finishes
#savefile outputfile1 outputfile2
for (( i=0; i<${#T_start1[@]}; i++))
do
    savefile ${outputs[i]}
done

########################################################
# Run the application, and we typically time it:
##------------------------------------------------------

echo "Handling control to julia script."
cd $SCRATCH
for (( i=0; i<${#T_start1[@]}; i++ ))
do
    srun -N1 --ntasks=1 $JULIA_PATH -p $CPUS $JULIA_SCRIPT ${T_start1[i]} ${T_start2[i]} ${T_start3[i]} ${T_init1[i]} ${T_init2[i]} ${T_init3[i]} > ${outputs[i]} &
done

wait
echo "Job script finished without issue."

#########################################################
# That was about all this time; lets call it a day...
##-------------------------------------------------------
# Finish the script
exit 0
