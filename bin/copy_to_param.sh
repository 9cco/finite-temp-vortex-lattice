#!/bin/bash
# Author: F.N. Krohg
# Created: 22.08.19
# Script for re-organizing simulation runs into forders of separate parameters.

gs=(0.3 0.5 1.0)
nus=(0.0 0.3 0.5 1.0)

for g in "${gs[@]}"
do
    for nu in "${nus[@]}"
    do
        for i in T\=*; do if [ -d "$i"/*g\=${g}_nu\=${nu}_fL\=1.0* ]; then cp -ur "$i"/*g\=${g}_nu\=${nu}_fL\=1.0* ../Parameter_Sort/g\=${g}/nu\=${nu}/; fi; done
    done
done
