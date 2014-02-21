#!/bin/bash

[[ -z $COMPILER ]] && COMPILER="gcc"

if [[ $COMPILER="gcc" ]]; then
   make -f Makefile_opteron
   echo "Running opteron tests"
   for size in {16384..49152..512}; 
   do
      ./dsymv_opteron $size
   done
   make -f Makefile_opteron clean
fi
