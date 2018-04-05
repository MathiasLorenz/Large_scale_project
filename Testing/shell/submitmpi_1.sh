#!/bin/sh 

# --  General options 

# Naming of the job and queue name
#BSUB -J mpi_1
#BSUB -q hpc

# Specify
#BSUB -oo Output.out 
#BSUB -eo Error.err 

# -- Technical options

# Ask for n cores placed on R host.
#BSUB -n 2
#BSUB -R "span[ptile=1]"

# Memory specifications. Amount we need and when to kill the
# program using too much memory.
#BSUB -R "rusage[mem=2GB]"
#BSUB -M 3GB

# Time specifications (hh:mm)
#BSUB -W 00:10

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
#echo LSB: job identifier is $LSB_JOBID
#echo LSB: executing queue is $LSB_QUEUE
#echo LSB: number of nodes is $LSB_NODE
#echo LSB: number of processors per node is $LSB_NUM_PPN
#echo LSB: total number of processors is $LSB_PROC
#echo LSB: working directory is $LSB_O_WORKDIR
#echo LSB: current home directory is $LSB_OUTDIR
echo LSB_OUTDIR: $LSB_OUTDIR
echo LSB_DJOB_NUMPROC: $LSB_DJOB_NUMPROC
echo LSB_MAX_NUM_PROCESSORS: $LSB_MAX_NUM_PROCESSORS
echo LSB_HOSTS: $LSB_HOSTS
echo LSB_MCPU_HOSTS: $LSB_MCPU_HOSTS
echo LSB_QUEUE: $LSB_QUEUE

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

	# Run the program
	N="10 20 50"
	for n in $N 
	do
		OUTPUT_INFO=full_matrix_mpi_z_slice \
			mpiexec -q -n $LSB_DJOB_NUMPROC ./jacobiSolver.bin mpi3d_1 $n \
			>> $LSB_JOBNAME-$n.dat

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
