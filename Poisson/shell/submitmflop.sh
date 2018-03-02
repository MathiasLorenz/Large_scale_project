#!/bin/sh

# ========== HOW TO ========== #
# To use this submition. Change the
# PBS name, the Meta variables and
# the modules needed ect. All output
# files from functions must be named
# outputSOMETHING.EXTENSION
# otherwise they are removed.

# Definition for the job and output
#PBS -N mflop
#PBS -o $PBS_JOBNAME.out
#PBS -e	$PBS_JOBNAME.err

# Definition of the technical stuff
#PBS -q hpcintro
#PBS -l walltime=1:00:00
#PBS -l feature=XeonE5-2660
#PBS -l nodes=1:ppn=1

# Meta variables
OUTPTH=../data
FIGPTH=../figures
MTLB=../matlabcode
LOGS=../logs

# Load needed modules
module load studio

# ======= GENERIC PART ======= #
# Directory and cpu info
cd $PBS_O_WORKDIR
lscpu >> $PBS_JOBNAME.cpu

# Run the actual code
echo Running $PBS_JOBNAME
./$PBS_JOBNAME.sh

# Export data with matlab
matlab -r "addpath('$OUTPTH');addpath('$MTLB');matlab$PBS_JOBNAME();ExportFigures();exit;"

# Cleanup
echo Cleaning
mv -f -t $LOGS $PBS_JOBNAME.out $PBS_JOBNAME.err $PBS_JOBNAME.cpu
mv -f -t $OUTPTH output*
mv -f -t $FIGPTH Figures/*
cd ../
rm -f -r $PBS_JOBNAME
