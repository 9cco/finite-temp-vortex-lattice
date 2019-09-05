#!/bin/bash
# Author: F.N. Krohg
# Created: 03.09.19
# Script for converting all the progress shots of S+ into movies.

Ts=(0.4 0.6 0.8 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 3.0)
progress_path="S+_progress"

for T in "${Ts[@]}"
do
    cd "T=${T}"
    for i in full_model_L\=32_*
    do
        cd "$i"
        if [ -d "$progress_path/T=$T" ]
        then
            cd "$progress_path/T=$T"
            if [ -f "step=0001.png" ]
            then
                ffmpeg -start_number 1 -framerate 24 -i step=%04d.png -pix_fmt yuv420p -c:v libx264 out.mp4
                rm *.png
            fi
            cd ../..
        else
            echo "No progress folder found in `pwd`"
        fi
        cd ..
    done
    cd ..
done
