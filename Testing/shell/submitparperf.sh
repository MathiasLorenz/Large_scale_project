#!/bin/sh 

# --  General options 

# Naming of the job and queue name
#BSUB -J parperf
#BSUB -q gpuv100

# Specify
#BSUB -oo Output.out 
#BSUB -eo Error.err 

# -- Technical options

# Ask for n cores placed on R host.
#BSUB -n 10
#BSUB -R "span[ptile=2]"

# Memory specifications. Amount we need and when to kill the
# program using too much memory.
#BSUB -R "rusage[mem=20GB]"
#BSUB -M 20GB

# Time specifications (hh:mm)
#BSUB -W 03:00

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
	module load cuda/9.1 mpi/2.1.0-gcc-6.3.0

}

# End of Preparation
#==============================================================================
# Define Program

Program()
{
	echo ' '
	echo Running computations

	start=`date +%s`
	# -------------------------------------------------------------------------
	# Define the actual test part of the script 

	#N="32"
	#C="2 3 4"

	C="2 3 4 5 6 7 8 9 10 11 12"
	N="512"

	# Run the MPI based tests
	TEST="mixed_1 mixed_2 mixed_3"
	for t in $TEST
	do
		for c in $C
		do
			dat="$c-$t.dat"
			echo "Test: $t, C: $c"
			mpiexec -q -n $c ./jacobiSolver.bin $t $N >> $LSB_JOBNAME-$dat
		done
	done

	# -------------------------------------------------------------------------
	end=`date +%s`

	runtime=$((end-start))
	echo "Time spent on computations: $runtime"
	mv -t $DPATH *.dat 
}

# End of Program
#==============================================================================
# Define Visualize

Visualize()
{
	echo ' '
	echo Visualizing
	matlab -r "addpath(genpath('../../'));matlabparperf('$LSB_JOBNAME','$DPATH/','$FIGS');exit;"
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
