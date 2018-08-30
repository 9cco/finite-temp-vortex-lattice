#!/bin/bash
# Author: F.N. Krohg
# Last updated: 30.08.18

error_exit()
{
	echo "$1" 1>&2
	exit 1
}

GET_SCRIPT="get_last_data.sh"
WRITE_SCRIPT="write_latex.sh"

# By recursively going through all symlinks of how the script was called
# we find the location of the script.
SCRIPT=$(readlink -f "$0")
# Then we use dirname to find the directory of this location.
SCRIPT_DIR=$(dirname "$SCRIPT")

cd $SCRIPT_DIR

[ -e $GET_SCRIPT ] || error_exit "ERROR: Could not find $GET_SCRIPT"
[ -e $WRITE_SCRIPT ] || error_exit "ERROR: Could not find $WRITE_SCRIPT"

echo "Getting latest data."
./$GET_SCRIPT
echo "Compiling LaTeX output."
./$WRITE_SCRIPT
