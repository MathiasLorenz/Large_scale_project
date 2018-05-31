#!/bin/sh 

# --  General options 

# Naming of the job and queue name
#BSUB -J profile
#BSUB -q gpum2050

# Specify
#BSUB -oo Output.out 
#BSUB -eo Error.err 

# -- Technical options

# Ask for n cores placed on R host.
#BSUB -n 2
#BSUB -R "span[ptile=1]"

# Memory specifications. Amount we need and when to kill the
# program using too much memory.
#BSUB -R "rusage[mem=20GB]"
#BSUB -M 20GB

# Time specifications (hh:mm)
#BSUB -W 00:10

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

	/appl/cuda/9.1/samples/bin/x86_64/linux/release/deviceQuery >> $LSB_JOBNAME.gpu
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

	# Run the programs (Max array size for GPU: 874)

	mpiexec -q -n $LSB_DJOB_NUMPROC nvprof \
		--output-profile profile.%q{OMPI_COMM_WORLD_RANK} \
		--process-name "rank %q{OMPI_COMM_WORLD_RANK}" \
		--context-name "rank %q{OMPI_COMM_WORLD_RANK}" \
		./jacobiSolver.bin mixed_2 100

	# -------------------------------------------------------------------------
	end=`date +%s`

	runtime=$((end-start))
	echo "Time spent on computations: $runtime"
	#mv -t $DPATH *.dat 
}

# End of Program
#==============================================================================
# Define Visualize

Visualize()
{
	echo ' '
	echo Visualizing
	#matlab -r "addpath(genpath('../../'));matlabperformance('$LSB_JOBNAME','$DPATH/','$FIGS');exit;"
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
