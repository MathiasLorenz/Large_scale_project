#!/bin/sh

# Define all needed folders relative to project head.
EPATH=Poisson
DPATH=Testing/data
FPATH=Testing/figures
LPATH=Testing/logs
SPATH=Testing/shell

# Make sure the excecutable is up to date
module load cuda/9.1 mpi/2.1.0-gcc-6.3.0
cd $EPATH; make realclean; make -s; 
if [ -f jacobiSolver.bin ]
then 
	cd ../
else
	echo "ERROR:"
	echo "	jacobiSolver.bin not found. Aborting tests."
	exit
fi

# Define files needed by the execution in all tests
EXEC="$EPATH/jacobiSolver.bin"

# Define which shell script will be executed
if [ -z "$1" ] ; then
	echo "ERROR:"
	echo "	Please specify which test to run."
	exit
else 
	TEST="$@"
fi

# Create all needed folders
mkdir -p $DPATH $FPATH $LPATH
echo Submitting the following tests:
echo ' '
echo $TEST
echo ' '

for test in $TEST
do
	# Create the folder needed
	rm -fr $LPATH/$test/*
	mkdir -p $LPATH/$test

	# Clean and copy all files needed
	cp -ft $LPATH/$test $SPATH/submit$test.sh $EXEC

	# Move to the directory submit the code and return
	cd $LPATH/$test
	bsub < submit$test.sh
	cd ../../../
done
