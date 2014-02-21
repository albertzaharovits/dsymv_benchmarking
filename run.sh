#!/bin/bash
#
# Author: Heri
#
# Script de submitere a job-urilor pe fiecare coda, folosind compilatoare diferite
#
#--modules "compilers/gcc-4.4.6" \

mprun.sh --job-name DSYMV-NEHALEM --queue ibm-nehalem.q \
	--script exec_script_nehalem.sh --show-qsub --show-script --batch-job
mprun.sh --job-name DSYMV-OPTERON --queue ibm-opteron.q \
	--script exec_script_opteron.sh --show-qsub --show-script --batch-job
mprun.sh --job-name DSYMV-QUAD --queue ibm-quad.q \
	--script exec_script_quad.sh --show-qsub --show-script --batch-job

