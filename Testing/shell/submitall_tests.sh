#!/bin/sh 

# --  General options 

# Naming of the job and queue name
#BSUB -J all_tests
#BSUB -q gpuv100

# Specify
#BSUB -oo Output.out 
#BSUB -eo Error.err 

# -- Technical options

# Add progress
#BSUB -K

# Ask for n cores placed on R host.
#BSUB -n 2
#BSUB -R "span[ptile=2]"

# Memory specifications. Amount we need and when to kill the
# program using too much memory.
#BSUB -R "rusage[mem=10GB]"
#BSUB -M 10GB

# Time specifications (hh:mm)
#BSUB -W 01:00

# GPU options
#BSUB -gpu "num=2:mode=exclusive_process"

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
	module load cuda/9.2 mpi/3.1.1-gcc-6.4.0-cuda-9.2
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

	export USE_TOLERANCE=on
	export OUTPUT_INFO=error

	# Run the program
	N="32 128"
	#N="10 20 50 100 200 300 400 500 600 700"
	for n in $N 
	do
		echo omp2d
		time OMP_NUM_THREADS=$LSB_DJOB_NUMPROC ./jacobiSolver.bin omp2d $n
		
		echo omp3d
		time OMP_NUM_THREADS=$LSB_DJOB_NUMPROC ./jacobiSolver.bin omp3d $n
		
		echo mpi3d_1
		time mpiexec -q -n 2 ./jacobiSolver.bin mpi3d_1 $n
		
		echo mpi3d_2
		time mpiexec -q -n $LSB_DJOB_NUMPROC ./jacobiSolver.bin mpi3d_2 $n
		
		echo mpi3d_3
		time mpiexec -q -n $LSB_DJOB_NUMPROC ./jacobiSolver.bin mpi3d_3 $n
		
		echo cuda_1
		time ./jacobiSolver.bin cuda_1 $n

		echo cuda_2
		time ./jacobiSolver.bin cuda_2 $n
		
		echo mixed_1
		time mpiexec -q -n $LSB_DJOB_NUMPROC ./jacobiSolver.bin mixed_1 $n
		
		echo mixed_2
		time mpiexec -q -n $LSB_DJOB_NUMPROC ./jacobiSolver.bin mixed_2 $n
		
		echo mixed_3
		time mpiexec -q -n $LSB_DJOB_NUMPROC ./jacobiSolver.bin mixed_3 $n
	
		echo mixed_4
		time mpiexec -q -n $LSB_DJOB_NUMPROC ./jacobiSolver.bin mixed_4 $n
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
	#matlab -r "addpath(genpath('../../'));matlab3Dplots('$LSB_JOBNAME','$DPATH/','$FIGS');exit;"
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
