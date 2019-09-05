#!/bin/bash
# Author: F.N. Krohg
# Last-updated: 27.05.19
# Pushes the content of neiboring Source, Scripts and bin directory to the project directory on Fram

error_exit()
{
	echo "$1" 1>&2
	exit 1
}


HOSTNAME="fram.sigma2.no"
UNAME="fredrkro"
PROJECT="nn2819k"
HOST_PATH_SOURCE="/cluster/projects/${PROJECT}/finite-temp-vortex-lattice/Source"
HOST_PATH_SCRIPTS="/cluster/projects/${PROJECT}/finite-temp-vortex-lattice/Scripts"
HOST_PATH_BIN="/cluster/projects/$PROJECT/finite-temp-vortex-lattice/bin"

# Make sure that the directory structure exists on host.
ssh $UNAME@$HOSTNAME "[ -d $HOST_PATH_SOURCE ] || mkdir -p $HOST_PATH_SOURCE; [ -d $HOST_PATH_SCRIPTS ] || mkdir -p $HOST_PATH_SCRIPTS; [ -d $HOST_PATH_BIN ] || mkdir -p $HOST_PATH_BIN"

# Finding source directory which should be neighboring to the script directory
# By recursively going through all symlinks of how the script was called
# we find the location of the script.
SCRIPT=$(readlink -f "$0")
# Then we use dirname to find the directory of this location.
BIN_DIR=$(dirname "$SCRIPT")

SOURCE_DIR="${BIN_DIR}/../Source"
SCRIPTS_DIR="${BIN_DIR}/../Scripts"
# Check that we have found Source directory
[ -d $SOURCE_DIR ] || error_exit "ERROR: Could not find Source directory in $SOURCE_DIR"
# Check that we have found bin directory
[ -d $BIN_DIR ] || error_exit "ERROR: Could not find bin directory in $BIN_DIR"
# Check that we have found Scripts directory
[ -d $SCRIPTS_DIR ] || error_exit "ERROR: Could not find bin directory in $SCRIPTS_DIR"

# Copy source files from Source dir to remote Source dir

echo "Copying $SOURCE_DIR to $HOSTNAME"
rsync -avz -e "ssh" --progress $SOURCE_DIR/ $UNAME@$HOSTNAME:$HOST_PATH_SOURCE
#echo "Copying $BIN_DIR to $HOSTNAME"
#rsync -avz -e "ssh" --progress $BIN_DIR/ $UNAME@$HOSTNAME:$HOST_PATH_BIN
echo "Copying $SCRIPTS_DIR to $HOSTNAME"
rsync -avz -e "ssh" --progress $SCRIPTS_DIR/ $UNAME@$HOSTNAME:$HOST_PATH_SCRIPTS

exit 0
