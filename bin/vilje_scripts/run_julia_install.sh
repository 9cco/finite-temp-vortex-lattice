#!/bin/bash
# Author: F.N. Krohg
# Last-updated: 30.08.18
# PBS job script for launching chiral superconductivity calculation on the vilje
# supercomputer. Based on a sample job script by Jabir Ali Ouassou.

##################################################################
# Configure the PBS queue system
##################################################################

#PBS -M fredrik.n.krohg@ntnu.no
#PBS -N install
#PBS -A nn2819k
#PBS -q workq
#PBS -o install.o
#PBS -e install.e
#PBS -l select=1:ncpus=1
#PBS -l walltime=00:15:00
#PBS -l pmem=1000MB

##################################################################
# Prepare the simulation
##################################################################

SOURCE_PATH="/work/fredrkro/finite-temp-vortex-lattice/Source"
JULIA_PATH="/home/ntnu/fredrkro/bin/julia"
JULIA_SCRIPT="install_packages.jl"
CPUS="1" # Needs to be the same as ncpus above

error_exit()
{
	echo "$1" 1>&2
	exit 1
}

# TODO: Use symbolic links for navigation references.

# Making sure we are in the right directory
[ -d $SOURCE_PATH ] || error_exit "ERROR: Could not find Source directory in $SOURCE_PATH"
cd $SOURCE_PATH
[ -f ./$JULIA_SCRIPT ] || error_exit "ERROR: Could not find julia script $JULIA_SCRIPT in `pwd`"

echo "Handling control to julia install script."
$JULIA_PATH $JULIA_SCRIPT
echo "Job script finished without issue."
