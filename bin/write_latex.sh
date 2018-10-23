#!/bin/bash
# Author: F.N. Krohg
# Last-updated: 30.08.18
# Assumes the current folder is full of files called somethign_plot.pdf which should be compiled
# into a single pdf file as well as a file called system_values.data containg information from
# the Monte-Carlo simulation. Takes these files and compiles a nice looking output pdf-file
# which is then opened in evince.

error_exit()
{
	echo "$1" 1>&2
	exit 1
}

# Read the simulation values from a file.

SYSTEM_FILE="system_values.data"
TEX_FILE="test_output.tex"
FOLDER_FILE="last_folder.txt"
LATEX_LOG="latex_compile.log"
JULIA_OUTPUT="julia.output"


# Figuring out if a command line argument was given
if [ -z $1 ]
then
	# If no command line arguments were given

	# By recursively going through all symlinks of how the script was called
	# we find the location of the script.
	SCRIPT=$(readlink -f "$0")
	# Then we use dirname to find the directory of this location.
	SCRIPT_DIR=$(dirname "$SCRIPT")
	# Position ourselves in the neighboring Data directory

	[ -d "$SCRIPT_DIR/../Data" ] || error_exit "ERROR: Could not find directory $SCRIPT_DIR/../Data"
	cd $SCRIPT_DIR/../Data
	# Check that last_folder.txt is in the Data directory.
	[ -e "./$FOLDER_FILE" ] || error_exit "ERROR: Could not find $FOLDER_FILE"
	LAST_FOLDER=`cat $FOLDER_FILE`
	[ -d "./$LAST_FOLDER" ] || error_exit "ERROR: `pwd`/$LAST_FOLDER is not a directory"
	cd $LAST_FOLDER
else
	# Check is argument exists.
	[ -d "./$1" ] || error_exit "ERROR: `pwd`$1 is not a directory"
	cd $1
fi

# Now we are in the directory where we assume the pdf-files and system_values.data is

# Try to find any file that ends on julia.output
julia_outputs=(`ls -1 | grep ".*$JULIA_OUTPUT"`)
[ ${#julia_outputs[@]} == 1 ] || error_exit "ERROR: Too many files ending with *$JULIA_OUTPUT"
JULIA_OUTPUT="${julia_outputs[0]}"

# Checking if tex_file exists and if it does, ask if we can overwrite
if [ -e $TEX_FILE ]
then
	echo -n "File $TEX_FILE already exists. Ok to overwrite? [y/n]: "
	read answ
	if [ $answ = "n" ]
	then
		exit 1
	fi
fi

[ -e "$SYSTEM_FILE" ] || error_exit "ERROR: Could not find $SYSTEM_FILE in `pwd`"
echo -e "Reading system constants.\n"
i=0
declare -a keywords
declare -a values
while read -r line
do
	keywords[i]=`expr "$line" : '^\([a-zA-Z_]*\)[\ ]*.*$'`
	values[i]=`expr "$line" : '^[a-zA-Z_]*[\ ]*\(.*\)$'`

	echo "${keywords[$i]} = ${values[$i]}"

	case ${keywords[i]} in
		"L")
			L=${values[i]}
			;;
		"GAMMA")
			GAMMA=${values[i]}
			;;
		"g")
			G=${values[i]}
			;;
		"NU")
			NU=${values[i]}
			;;
		"f")
			F=${values[i]}
			;;
		"TEMP")
			TEMP=${values[i]}
			;;
		"INV_TEMP")
			INV_TEMP=${values[i]}
			;;
		"NR_MEASUREMENTS")
			NR_MEASUREMENTS=${values[i]}
			;;
		"MEASUREMENT_INTERVAL")
			MEASUREMENT_INTERVAL=${values[i]}
			;;
		"SIM_THETA_MAX")
			SIM_THETA_MAX=${values[i]}
			;;
		"SIM_UMAX")
			SIM_UMAX=${values[i]}
			;;
		"SIM_AMAX")
			SIM_AMAX=${values[i]}
			;;
        "THERM_T")
            THERM_T=${values[i]}
            ;;
		*)
			echo "${keywords[i]} not defined as variable."
			;;
	esac

	i=$(( $i + 1 ))
done < "$SYSTEM_FILE"

# Check that julia output-file is included.
[ -e $JULIA_OUTPUT ] || error_exit "ERROR: Could not find julia output in `pwd`/$JULIA_OUTPUT"

# Write the beginning of the tex file with correct values.
echo -e "\nWriting these constants to LaTeX file."

cat << EOF > $TEX_FILE || error_exit "ERROR: Can't write LaTeX file"
\documentclass[a4paper, 11pt]{article}

\usepackage{float}
\usepackage{graphicx}
%\usepackage[utf8]{inputenc}
%\usepackage[margin=2.5cm, bottom=0.75in, a4paper]{geometry}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{listingsutf8}
\usepackage{fancyvrb}

\usepackage{fontspec}
\newfontfamily\codefont{Unifont}

%\input{../../lib/functions.tex}

\title{Output of Markov-Chain Monte-Carlo simulation}
\author{Krohg, F.N.}
\date{\today}       

\begin{document}
\maketitle
\section{Program input}
\begin{table}[h]
  \centering
  \caption{Input parameters used in the simulation. \$L$ is the system size in one dimension.}
  \begin{tabular}{l l}\hline
	\$L$ & $L\\\\
	\$\gamma$ & $GAMMA\\\\
	\$g$ & $G\\\\
	\$\nu$ & $NU\\\\
	\$f$ & $F\\\\
	\$T$ & $TEMP\\\\
	\$\beta$ & $INV_TEMP\\\\
	\$M$ & $NR_MEASUREMENTS\\\\
    \$t_0$ & $THERM_T\\\\
	\$\Delta t$ & $MEASUREMENT_INTERVAL\\\\
	\$\theta_\text{max}$ & $SIM_THETA_MAX\\\\
	\$u_\text{max}$ & $SIM_UMAX\\\\
	\$A_\text{max}$ & $SIM_AMAX\\\\\hline
  \end{tabular}
  \label{tab:input_parameters}
\end{table}
\clearpage

\section{Output}
\VerbatimInput[frame=lines, label=Julia Output, fontfamily=Unifont(0), fontsize=\scriptsize, framesep=5mm]{$JULIA_OUTPUT}
%\lstinputlisting{$JULIA_OUTPUT}
\clearpage

\section{Plots}
EOF

# Find all the filenames ending with *_plot.pdf and include then in latex as a figure.

# Getting list of filenames we want to include
# We use ls to ceate a list of all files in the current directory, then pipe this to grep
# which filters the list to only contain files ending in _plot.pdf using a regular expression
# and finally convert the space-delimitered string of files returned from grep into an array
# in bash. Not that this failes if one of the filenames ending in _plot.pdf contains a space,
# since the way bash defines arrays is space-deilimitered.
include_files=(`ls -1 | grep ".*_plot\.pdf"`)

# Loop over the files in include_files
echo "Including *_plots.pdf files as figures in LaTeX document."
for file in "${include_files[@]}"; do
	# Strip the _plot.pdf from the filename using a regular expression.
	echo "Adding $file"
	LABEL="`expr $file : '\(.*\)_plot.pdf'`"
	# Write the included file into a figure
	cat << EOF >> $TEX_FILE || error_exit "ERROR: Can't write LaTeX file"
\begin{figure}[h]
  \centering
  \includegraphics[width=\textwidth]{$file}
  \label{fig:$LABEL}
\end{figure}
EOF
done

# Write the final part of the latex document. The -e is added to echo to turn of interpretation
# of backslashes so that we can write newline characters.
echo -e "\n\\\end{document}\n" >> $TEX_FILE

echo "Compiling LaTeX document."
# The we compile the latex document
xelatex $TEX_FILE > $LATEX_LOG 2>&1

# And open the created pdf document. We do this by taking the .tex file-name and stripping the
# trailing .tex with .pdf using a regular expression. The pdf is opened as a background process
# and all output is discarded.
evince "`expr $TEX_FILE : '\(.*\)\.tex'`.pdf" 1>/dev/null 2>&1 &
