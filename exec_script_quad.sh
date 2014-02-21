#!/bin/bash

[[ -z $COMPILER ]] && COMPILER="gcc"

if [[ $COMPILER="gcc" ]]; then
   make -f Makefile_quad
   echo "Running quad tests"
   for size in {16384..49152..512}; 
   do
      ./dsymv_quad $size
   done
   make -f Makefile_quad clean
fi
