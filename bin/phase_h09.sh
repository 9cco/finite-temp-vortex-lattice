#!/bin/bash
# PBS job script for launching superconductivity simulations. 
# Written by Jabir Ali Ouassou <jabirali@switzerlandmail.ch>
# Based on an example script by Roy Dragseth and Egil Holsvik.

##################################################################
# Configure the PBS queue system
##################################################################
#PBS -M jabir.a.ouassou@ntnu.no
#PBS -A ntnu652
#PBS -q workq
#PBS -lselect=1:ncpus=32
#PBS -lwalltime=48:00:00
#PBS -lpmem=1000MB

##################################################################
# Prepare the simulation
##################################################################

# Specify the simulation title and parameters
title="ZSV_h09"

# Specify the parameters used to launch programs
parameters=( 
  '0.9 0.00 0.01'
  '0.9 0.00 1.01'
  '0.9 0.05 0.01'
  '0.9 0.05 1.01'
  '0.9 0.10 0.01'
  '0.9 0.10 1.01'
  '0.9 0.15 0.01'
  '0.9 0.15 1.01'
  '0.9 0.20 0.01'
  '0.9 0.20 1.01'
  '0.9 0.25 0.01'
  '0.9 0.25 1.01'
  '0.9 0.30 0.01'
  '0.9 0.30 1.01'
  '0.9 0.35 0.01'
  '0.9 0.35 1.01'
  '0.9 0.40 0.01'
  '0.9 0.40 1.01'
  '0.9 0.45 0.01'
  '0.9 0.45 1.01'
  '0.9 0.50 0.01'
  '0.9 0.50 1.01'
  '0.9 0.55 0.01'
  '0.9 0.55 1.01'
  '0.9 0.60 0.01'
  '0.9 0.60 1.01'
  '0.9 0.65 0.01'
  '0.9 0.65 1.01'
  '0.9 0.70 0.01'
  '0.9 0.70 1.01'
  '0.9 0.75 0.01'
  '0.9 0.75 1.01'
  '0.9 0.80 0.01'
  '0.9 0.80 1.01'
  '0.9 0.85 0.01'
  '0.9 0.85 1.01'
  '0.9 0.90 0.01'
  '0.9 0.90 1.01'
  '0.9 0.95 0.01'
  '0.9 0.95 1.01'
  '0.9 1.00 0.01'
  '0.9 1.00 1.01'
  )

# Load any required modules
module load intelcomp/18.0.0

# Add the required programs to the system path
PATH="~/Code/bin:$PATH"


##################################################################
# Perform the simulation
##################################################################

# Launch the simulation script once for each parameter defined above
for n in $(seq 0 $(( ${#parameters[@]} - 1))); do
  # Make a work directory
  workdir="/work/$LOGNAME/$(date +'%Y-%m-')${title}/${n}"
  mkdir -p "$workdir"
  cd "$workdir"

  # Execute the simulation script in the background
  eval converge $PBS_O_WORKDIR/materials.conf ${parameters[n]} &
done

# Wait for all tasks to finish
echo "Waiting for tasks to complete"
wait
echo "All jobs complete!"
