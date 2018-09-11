#!/bin/bash
# Author: F.N. Krohg
# Last-updated: 30.08.18
# PBS job script for launching chiral superconductivity calculation on the vilje
# supercomputer. Based on a sample job script by Jabir Ali Ouassou.

##################################################################
# Configure the PBS queue system
##################################################################

#PBS -M fredrik.n.krohg@ntnu.no
#PBS -N sfvl-test
#PBS -A nn2819k
#PBS -q workq
#PBS -o vortex.o
#PBS -e vortex.e
#PBS -l select=1:ncpus=16
#PBS -l walltime=01:00:00
#PBS -l pmem=1000MB

##################################################################
# Prepare the simulation
##################################################################

SOURCE_PATH="/work/fredrkro/finite-temp-vortex-lattice/Source"
JULIA_PATH="/home/ntnu/fredrkro/bin/julia"
JULIA_SCRIPT="driver.jl"
CPUS="15" # Needs to be the same as ncpus above
OUTPUT="julia.output"
FOLDER_FILE="last_folder.txt"
REL_DATA_PATH="../Data"

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

# Setup variables
g="0.3"
NU="0.3"
H="-0.72"
L="100"
TEMP="0.05"
GAMMA="1.0"
M="1000"
dt="3000"

echo "Handling control to julia script. Writing output to $OUTPUT"
$JULIA_PATH -p $CPUS $JULIA_SCRIPT $g $NU $H $L $TEMP $GAMMA $M $dt > $OUTPUT 2>&1

# Finding output directory
[ -d $REL_DATA_PATH ] || error_exit "ERROR: Could not find neighboring Data directory."
[ -e "$REL_DATA_PATH/$FOLDER_FILE" ] || error_exit "ERROR: Could not find $FOLDER_FILE in `pwd`"
LAST_FOLDER=`cat $REL_DATA_PATH/$FOLDER_FILE`
[ -d "$REL_DATA_PATH/$LAST_FOLDER" ] || error_exit "ERROR: `pwd`/$REL_DATA_PATH/$LAST_FOLDER is not a directory"

# Copying output file to data directory
echo "Copying output of $OUTPUT to $REL_DATA_PATH/$LAST_FOLDER/"
cp $OUTPUT $REL_DATA_PATH/$LAST_FOLDER/

echo "Job script finished without issue."
