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
#SBATCH --job-name=MCMCNo
## Allocating amount of resources:
#SBATCH --nodes=4
## Number of tasks (aka processes) to start on each node: Pure mpi, one task per core
#SBATCH --ntasks-per-node=32 --cpus-per-task=1
## No memory pr task since this option is turned off on Fram in QOS normal.
## Run for x minutes, syntax is d-hh:mm:ss
#SBATCH --time=2-00:00:00 

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

error_exit()
{
	echo "$1" 1>&2
	exit 1
}

SOURCE_PATH="/cluster/projects/nn2819k/finite-temp-vortex-lattice/Source/"
SCRIPTS_PATH="/cluster/projects/nn2819k/finite-temp-vortex-lattice/Scripts/"
JULIA_PATH="/cluster/home/fredrkro/Programs/julia-1.0.4/bin/julia"

# We use 4 nodes and will thus run 4 simulations for 4 different temperatures
# which should output to 4 diffferent output files
gs=(1.0 1.0 1.0 1.0 0.3 0.3 0.3 0.3 0.5 0.5 0.5 0.5)
nus=(0.0 0.3 0.5 1.0 0.0 0.3 0.5 1.0 0.0 0.3 0.5 1.0)

# Check that arrays are equal
[ ${#gs[@]} = ${#nus[@]} ] || error_exit "ERROR: input arrays of un-equal length"

declare -a outputs
for (( i=0; i<${#nus[@]}; i++))
do
    outputs[i]="g${gs[i]}_nu${nus[i]}.out"
done

# Needs to be edited when changing scripts
JULIA_SCRIPT="driver_new.jl"
CPUS="8" # Max is --ntasks-per-node * --nodes


# Making sure we are in the right directory
[ -d $SCRIPTS_PATH ] || error_exit "ERROR: Could not find Source directory in $SCRIPTS_PATH"
cd $SCRIPTS_PATH
[ -f ./$JULIA_SCRIPT ] || error_exit "ERROR: Could not find julia script $JULIA_SCRIPT in `pwd`"
cp $JULIA_SCRIPT $SCRATCH/

#######################################################
## Prepare jobs, moving input files and making sure 
# output is copied back and taken care of
##-----------------------------------------------------

# Prepare input files
#cp inputfiles $SCRATCH
#cd $SCRATCH

# Make sure output is copied back after job finishes
#savefile outputfile1 outputfile2
for (( i=0; i<${#nus[@]}; i++))
do
    savefile ${outputs[i]}
done

########################################################
# Run the application, and we typically time it:
##------------------------------------------------------

echo "Handling control to julia script."
cd $SCRATCH
for (( i=0; i<${#nus[@]}; i++ ))
do
    $JULIA_PATH -p $CPUS $JULIA_SCRIPT ${gs[i]} ${nus[i]} > ${outputs[i]} &
done

wait
echo "Job script finished without issue."

#########################################################
# That was about all this time; lets call it a day...
##-------------------------------------------------------
# Finish the script
exit 0
