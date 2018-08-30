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
DOMAIN="vilje.hpc.ntnu.no"
USER="fredrkro"

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
echo "Downloading last folder name to `pwd`"
scp $USER@$DOMAIN:$EXTERNAL_TREE/$FOLDER_FILE ./ || error_exit "ERROR: Could not download file $FOLDER_FILE"
LAST_FOLDER=`cat $FOLDER_FILE`

# Checking if this folder already exists.
if [ -d "./$LAST_FOLDER" ]
then
	echo -n "The folder $LAST_FOLDER already exists. Overwrite? [y/n]: "
	read answ
	if [ $answ = "n" ]
	then
		exit 1
	fi
fi

echo "Downloading the last folder $LAST_FOLDER to `pwd`"
# Download the last used folder recursively.
scp -r $USER@$DOMAIN:$EXTERNAL_TREE/$LAST_FOLDER ./ || error_exit "ERROR: Could not download simulation data from $DOMAIN:$EXTERNAL_TREE/$LAST_FOLDER"
