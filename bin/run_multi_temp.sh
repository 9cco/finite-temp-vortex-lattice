#!/bin/bash
# Author: F.N. Krohg
# Last-updated: 11.09.18
# PBS job script for launching chiral superconductivity calculation on the vilje
# supercomputer. Based on a sample job script by Jabir Ali Ouassou.
# Launches numerous of the run-single-job.sh job

error_exit()
{
	echo "$1" 1>&2
	exit 1
}


# First we need to input the temperature range
TEMPS=(0.1109 0.1107 0.1105 0.1103 0.1101)
# Then the other variables are set as before
g="0.3"
NU="0.3"
H="-0.72"
L="80"
GAMMA="1.0"
M="300"
dt="5000"

# For each of the values in TEMPS, we make a new name
declare -a IDS
declare -a names
t_len=${#TEMPS[@]}
for (( i=1; i<=${t_len}; i++ ));
do
	names+=("sfvl3de_$i")
done

# For each of the temps we create a separate temp_single_job.pbs script
# which can be sent to the PBS scheduling system.
for (( i=0; i<$t_len; i++ ));
do
	cat << EOF > ${names[$i]}.pbs || error_exit "ERROR: Can't write PBS script $i"
#!/bin/bash
# Author: F.N. Krohg
# Last-updated: `date`
# PBS job script for launching chiral superconductivity calculation on the vilje
# supercomputer. Based on a sample job script by Jabir Ali Ouassou.

##################################################################
# Configure the PBS queue system
##################################################################

#PBS -M fredrik.n.krohg@ntnu.no
#PBS -N ${names[$i]}
#PBS -A nn2819k
#PBS -q workq
#PBS -o ${names[$i]}.o
#PBS -e ${names[$i]}.e
#PBS -l select=1:ncpus=16
#PBS -l walltime=62:00:00
#PBS -l pmem=1000MB

##################################################################
# Prepare the simulation
##################################################################

SOURCE_PATH="/work/fredrkro/finite-temp-vortex-lattice/Source"
DATA_PATH="/work/fredrkro/finite-temp-vortex-lattice/Data"
JULIA_PATH="/home/ntnu/fredrkro/bin/julia"
JULIA_SCRIPT="driver.jl"
CPUS="15" # Needs to be the same as ncpus above
OUTPUT="${names[$i]}_julia.output"
FOLDER_FILE="last_folder.txt"
NEED_LATEX_FILE="need_latex.txt"
REL_DATA_PATH="../Data"

error_exit()
{
	echo "\$1" 1>&2
	exit 1
}

# TODO: Use symbolic links for navigation references.

# Making sure necessary paths and script exist.
[ -d \$SOURCE_PATH ] || error_exit "ERROR: Could not find Source directory in \$SOURCE_PATH"
[ -f \$SOURCE_PATH/\$JULIA_SCRIPT ] || error_exit "ERROR: Could not find julia script \$JULIA_SCRIPT in \$SOURCE_PATH"

# Setup variables
g="$g"
NU="$NU"
H="$H"
L="$L"
TEMP="${TEMPS[$i]}"
GAMMA="$GAMMA"
M="$M"
dt="$dt"

# Setup work directory
[ -d \$DATA_PATH ] || error_exit "ERROR: Could not find Data directory in \$DATA_PATH"
cd \$DATA_PATH
# Make work directory name
WORK_NAME="VORTEX_GAMMA_\${GAMMA}_g_\${g}_NU_\${NU}_H_\${H}_T_\${TEMP}_L_\${L}_M_\${M}"
# Make directory if it does not exist.
echo "Entering directory \$WORK_NAME"
[ -d \$WORK_NAME ] || mkdir \$WORK_NAME
# Write directory name to FOLDER_FILE
echo "\$WORK_NAME" > \$FOLDER_FILE
# Add name to NEED_LATEX_FILE
[ -e \$NEED_LATEX_FILE ] || touch \$NEED_LATEX_FILE
echo "\$WORK_NAME" >> \$NEED_LATEX_FILE
# Change to directory
cd \$WORK_NAME

# Now we are in the correct work directory.

echo "Handling control to julia script. Writing output to \$OUTPUT"
\$JULIA_PATH -p \$CPUS \$SOURCE_PATH/\$JULIA_SCRIPT \$g \$NU \$H \$L \$TEMP \$GAMMA \$M \$dt > \$OUTPUT 2>&1


echo "Job script finished without issue."
EOF
	# Make the script executable
	chmod a+x ${names[$i]}.pbs
	# Send the script to the PBS scheduling system
	IDS+=("`qsub ${names[$i]}.pbs`")
	# Check the status of sent job
	qstat ${IDS[$i]}
	echo "Job started: ${IDS[$i]}"
done
echo "All jobs started, exiting"
