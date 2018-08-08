#!/bin/sh 

# --  General options 

# Naming of the job and queue name
#BSUB -J cuda_1
#BSUB -q gpuv100

# Specify
#BSUB -oo Output.out 
#BSUB -eo Error.err 

# -- Technical options

# Ask for n cores placed on R host.
#BSUB -n 4
#BSUB -R "span[ptile=1]"
#BSUB -K

# Memory specifications. Amount we need and when to kill the
# program using too much memory.
#BSUB -R "rusage[mem=20GB]"
#BSUB -M 20GB

# Time specifications (hh:mm)
#BSUB -W 01:00

# GPU options
#BSUB -gpu "num=1:mode=exclusive_process"

# -- Notification options

# Set the email to recieve to and when to recieve it
##BSUB -u your_email_address
#BSUB -B		# Send notification at start
#BSUB -N 		# Send notification at completion

echo --------------------------------------------------------------------------
echo 'Job: '$LSB_JOBNAME', is running on '$LSB_DJOB_NUMPROC' cores.'
echo --------------------------------------------------------------------------
echo LSB: job identifier is $LSB_JOBID
echo LSB: execution queue is $LSB_QUEUE
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
	# nvidia-smi
}

# End of Preparation
#==============================================================================
# Define Program

Program()
{
	echo ' '
	echo Running computations
	echo ' '
	# -------------------------------------------------------------------------
	# Define the actual test part of the script 
	# Run the program

	export MAX_ITER=10894
	
	N="128"
	for n in $N 
	do
		echo CUDA 1: Clean
			USE_TOLERANCE=off ./jacobiSolver.bin cuda_1 $n

		export USE_TOLERANCE=on

		echo CUDA 1: Tolerance
			./jacobiSolver.bin cuda_1 $n
		
		echo MPI 3: Error
			OUTPUT_INFO=error \
			mpiexec -qn $LSB_DJOB_NUMPROC ./jacobiSolver.bin mpi3d_3 $n

		echo CUDA 1: Error
			OUTPUT_INFO=error \
			./jacobiSolver.bin cuda_1 $n

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
	# matlab -r "addpath(genpath('../../'));matlab3Dplots('$LSB_JOBNAME','$DPATH/','$FIGS');exit;"
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
