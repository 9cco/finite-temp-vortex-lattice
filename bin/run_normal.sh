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

SOURCE_PATH="/cluster/projects/nn2819k/finite-temp-vortex-lattice/Source/"
SCRIPTS_PATH="/cluster/projects/nn2819k/finite-temp-vortex-lattice/Scripts/"
JULIA_PATH="/cluster/home/fredrkro/Programs/julia-1.0.4/bin/julia"

# We use 4 nodes and will thus run 4 simulations for 4 different temperatures
# which should output to 4 diffferent output files
gs=(0.31622776601683794 0.31622776601683794 0.31622776601683794 0.31622776601683794)
nus=(-0.3, 0.3, 0.6, 0.9)
OUTPUT1="g${gs[0]}_nu${nus[0]}.out"
OUTPUT2="g${gs[1]}_nu${nus[1]}.out"
OUTPUT3="g${gs[2]}_nu${nus[2]}.out"
OUTPUT4="g${gs[3]}_nu${nus[3]}.out"

# Needs to be edited when changing scripts
JULIA_SCRIPT="driver_new.jl"
CPUS="31" # Max is --ntasks-per-node * --nodes

error_exit()
{
	echo "$1" 1>&2
	exit 1
}

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
savefile $OUTPUT1 $OUTPUT2 $OUTPUT3 $OUTPUT4

########################################################
# Run the application, and we typically time it:
##------------------------------------------------------

echo "Handling control to julia script."
cd $SCRATCH
$JULIA_PATH -p $CPUS $JULIA_SCRIPT ${gs[0]} ${nus[0]} > $OUTPUT1 &
$JULIA_PATH -p $CPUS $JULIA_SCRIPT ${gs[1]} ${nus[1]} > $OUTPUT2 &
$JULIA_PATH -p $CPUS $JULIA_SCRIPT ${gs[2]} ${nus[2]} > $OUTPUT3 &
$JULIA_PATH -p $CPUS $JULIA_SCRIPT ${gs[3]} ${nus[3]} > $OUTPUT4 &
wait
echo "Job script finished without issue."

#########################################################
# That was about all this time; lets call it a day...
##-------------------------------------------------------
# Finish the script
exit 0
