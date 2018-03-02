#!/bin/sh 

# --  General options 

# Naming of the job and queue name
#BSUB -J My_Application
#BSUB -q hpc

# Specify
#BSUB -oo Output_%J.out 
#BSUB -eo Error_%J.err 

# -- Technical options

# Ask for 4 cores placed on 1 host.
#BSUB -n 4 
#BSUB -R "span[hosts=1]"

# Memory specifications. Amount we need and when to kill the
# program using too much memory.
#BSUB -R "rusage[mem=2GB]"
#BSUB -M 3GB

# Time specifications (hh:mm)
#BSUB -W 24:00 

# -- Notification options

# Set the email to recieve to and when to recieve it
##BSUB -u your_email_address
#BSUB -B		# Send notification at start
#BSUB -N 		# Send notification at completion
echo %J 
exit;

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
