#!/bin/bash
# Author: F.N. Krohg
# Last-updated: 30.08.18

error_exit()
{
	echo "$1" 1>&2
	exit 1
}

FOLDER_FILE="last_folder.txt"
EXTERNAL_TREE="/work/fredrkro/finite-temp-vortex-lattice/Data"
HOSTNAME="vilje.hpc.ntnu.no"
UNAME="fredrkro"

# By recursively going through all symlinks of how the script was called
# we find the location of the script.
SCRIPT=$(readlink -f "$0")
# Then we use dirname to find the directory of this location.
SCRIPT_DIR=$(dirname "$SCRIPT")
# Position ourselves in the neighboring Data directory

[ -d "$SCRIPT_DIR/../Data" ] || error_exit "ERROR: Could not find directory $SCRIPT_DIR/../Data"
cd $SCRIPT_DIR/../Data

# Find the last data directory
# Download the file containing the last updated folder name
echo "Downloading Data directory to `pwd`"
rsync -avx -e "ssh" --progress $UNAME@$HOSTNAME:$EXTERNAL_TREE/ .

