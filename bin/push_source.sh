#!/bin/bash
# Author: F.N. Krohg
# Last-updated: 09.09.18
# Pushes the content of neiboring Source and bin directory to the work directory on vilje

error_exit()
{
	echo "$1" 1>&2
	exit 1
}


HOSTNAME="vilje.hpc.ntnu.no"
UNAME="fredrkro"
PROJECT="finite-temp-vortex-lattice"
HOST_PATH_SOURCE="/work/${UNAME}/${PROJECT}/Source"
HOST_PATH_DATA="/work/${UNAME}/${PROJECT}/Data"
HOST_PATH_BIN="/work/$UNAME/$PROJECT/bin"

# Make sure that the directory structure exists on host.
ssh $UNAME@$HOSTNAME "[ -d $HOST_PATH_SOURCE ] || mkdir -p $HOST_PATH_SOURCE; [ -d $HOST_PATH_DATA ] || mkdir -p $HOST_PATH_DATA; [ -d $HOST_PATH_BIN ] || mkdir -p $HOST_PATH_BIN"

# Finding source directory which should be neighboring to the script directory
# By recursively going through all symlinks of how the script was called
# we find the location of the script.
SCRIPT=$(readlink -f "$0")
# Then we use dirname to find the directory of this location.
SCRIPT_DIR=$(dirname "$SCRIPT")

SOURCE_DIR="${SCRIPT_DIR}/../Source"
# Check that we have found Source directory
[ -d $SOURCE_DIR ] || error_exit "ERROR: Could not find Source directory in $SOURCE_DIR"
# Check that we have found bin directory
[ -d $SCRIPT_DIR ] || error_exit "ERROR: Could not find bin directory in $SCRIPT_DIR"

# Copy source files from Source dir to remote Source dir

echo "Copying $SOURCE_DIR to $HOSTNAME"
rsync -avz -e "ssh" --progress $SOURCE_DIR/ $UNAME@$HOSTNAME:$HOST_PATH_SOURCE
echo "Copying $SCRIPT_DIR to $HOSTNAME"
rsync -avz -e "ssh" $SCRIPT_DIR/ $UNAME@$HOSTNAME:$HOST_PATH_BIN

exit 0
