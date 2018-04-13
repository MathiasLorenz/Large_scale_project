#!/bin/sh 

# --  General options 

# Naming of the job and queue name
#BSUB -J mflop
#BSUB -q gpuv100

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
#BSUB -W 00:02 

# -- Notification options

# Set the email to recieve to and when to recieve it
##BSUB -u your_email_address
#BSUB -B		# Send notification at start
#BSUB -N 		# Send notification at completion

LSB_NODES=$(echo | grep -c '' $LSB_HOSTS)
LSB_PROC=$(($LSB_MAX_NUM_PROCESSORS/$LSB_NODES))
echo --------------------------------------------------------------------------
echo 'Job: '$LSB_JOBNAME', is running on '$LSB_NODES' nodes'
echo --------------------------------------------------------------------------
echo LSB: job identifier is $LSB_JOBID
echo LSB: execution queue is $LSB_QUEUE
echo LSB: number of nodes is $LSB_NODE
echo LSB: number of processors per node is $LSB_PROC
echo LSB: total number of processors is $LSB_MAX_NUM_PROCESSORS
echo LSB: working directory is $LSB_OUTDIR
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
	
	# -------------------------------------------------------------------------
	# Define the actual test part of the script 

	# Define program variables
	#N="15 20 25 30 40 50 100 200 400 800 1000 1500 2000 4000"
	N="10 50 100"

	MET="omp2d omp3d"

	# Run the program
	OUTPUT_INFO=timing
	USE_TOLERANCE=off
	MAX_ITER=1000

	for met in $MET
	do
		echo $met 
		for n in $N
		do
			echo $n
			./jacobiSolver.bin $met $n $n $n >> $LSB_JOBNAME$met.dat
		done
	done

	# -------------------------------------------------------------------------
	mv -t $DPATH *.dat 
}

# End of Program
#==============================================================================
# Define Visualize

Visualize()
{
	echo ' '
	echo Visualizing
	matlab -r "addpath(genpath('../../'));matlab$LSB_JOBNAME('$DPATH/','$FIGS');exit;"
}

# End of Visualize
#==============================================================================
# Define Finalize

Finalize()
{
	echo ' '
	echo Finalizing

	mv -ft $FPATH $FIGS/*

	echo Figures moved to $FPATH.
	echo Test concluded successfully.
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

echo ' '
echo --------------------------------------------------------------------------

Program

echo ' '
echo --------------------------------------------------------------------------

Visualize

echo ' '
echo --------------------------------------------------------------------------

Finalize

echo ' '
echo --------------------------------------------------------------------------
# ==============================   End of File   ==============================
