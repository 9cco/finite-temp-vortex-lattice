#!/bin/bash
# Author: F.N. Krohg
# Created: 22.08.19
# Script for running the post_temp_series.jl script within
# all the parameter folders.

gs=(0.3 0.5 1.0)
nus=(0.0 0.3 0.5 1.0)

for g in "${gs[@]}"
do
    cd "g=$g";
    for nu in "${nus[@]}"
    do
        cd "nu=$nu";
        julia -p 1 ~/mc/Scripts/post_temp_series.jl;
        cd ..;
    done;
    cd ..;
done
