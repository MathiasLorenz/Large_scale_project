#!/bin/sh 

# --  General options 

# Naming of the job and queue name
#BSUB -J mflop
#BSUB -q hpc

# Specify
#BSUB -oo Output.out 
#BSUB -eo Error.err 

# -- Technical options

# Ask for 4 cores placed on 1 host.
#BSUB -n 1 
#BSUB -R "span[hosts=1]"

# Memory specifications. Amount we need and when to kill the
# program using too much memory.
#BSUB -R "rusage[mem=2GB]"
#BSUB -M 3GB

# Time specifications (hh:mm)
#BSUB -W 00:01 

# -- Notification options

# Set the email to recieve to and when to recieve it
##BSUB -u your_email_address
#BSUB -B		# Send notification at start
#BSUB -N 		# Send notification at completion

LSB_PROC=$(echo | grep -c '' $LSB_NODEFILE)
# LSB_NODE=$(($LSB_PROC/$LSB_NUM_PPN))
echo --------------------------------------------------------------------------
echo 'Job: '$LSB_JOBNAME', is running on '$LSB_NODE' nodes'
echo --------------------------------------------------------------------------
echo LSB: job identifier is $LSB_JOBID
echo LSB: executing queue is $LSB_QUEUE
echo LSB: number of nodes is $LSB_NODE
echo LSB: number of processors per node is $LSB_NUM_PPN
echo LSB: total number of processors is $LSB_PROC
echo LSB: working directory is $LSB_O_WORKDIR
echo LSB: current home directory is $LSB_O_HOME
echo --------------------------------------------------------------------------

# End of LSB info
#==============================================================================
# Define Preparation

DPATH=../../data
FPATH=../../figures
LPATH=../
MPATH=../../matlab

LOGS="Output.out Error.err $LSB_JOBNAME.cpu"
DATA=$LSB_JOBNAME.dat
FIGS=$LSB_JOBNAME.fig

Prepare()
{
	echo ' '
	echo Preparing
	mkdir -p $FIGS
	lscpu >> $LSB_JOBNAME.cpu
	
	# Define modules
	module load cuda/9.1 mpi/2.1.0-gcc-6.3.0
}

# End of Preparation
#==============================================================================
# Define Program

Program()
{
	echo ' '
	echo Running computations
	./$LSB_JOBNAME.sh

	mv -t $DPATH *.dat 
}

# End of Program
#==============================================================================
# Define Visualize

Visualize()
{
	echo ' '
	echo Visualizing
	cd $MPATH
	matlab -r "matlab$LSB_JOBNAME('$FIGS');exit;"
	cd ../logs/$LSB_JOBNAME
}

# End of Visualize
#==============================================================================
# Define Finalize

Finalize()
{
	echo ' '
	echo Finalizing

	mv -ft $FPATH $FIGS/*
	
}

# End of Visualize
#==============================================================================
# Define Early Termination

early()
{
	echo ' '
	echo ' ================= WARNING: EARLY TERMINATION ================= '
	echo ' '
}
trap 'early' 2 9 15

# End of Early Termination
#==============================================================================
# Call Functions

Prepare
Program
Visualize
Finalize

echo ' '
echo --------------------------------------------------------------------------
# ==============================   End of File   ==============================
