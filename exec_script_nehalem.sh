#!/bin/bash

[[ -z $COMPILER ]] && COMPILER="gcc"

if [[ $COMPILER="gcc" ]]; then
   make -f Makefile_nehalem
   echo "Running nehalem tests"
   for size in {16384..49152..512}; 
   do
      ./dsymv_nehalem $size
   done
   make -f Makefile_nehalem clean
fi
