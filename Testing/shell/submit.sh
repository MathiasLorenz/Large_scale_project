#!/bin/sh

# Define all needed folders relative to project head.
EPATH=Poisson
DPATH=Testing/data
FPATH=Testing/figures
LPATH=Testing/logs
SPATH=Testing/shell

# Make sure the excecutable is up to date
cd $EPATH; make; cd ../

# Define files needed by the execution in all tests
FILE="$EPATH/jacobiSolver.bin"

# Define which shell script will be executed
VER="mflop"

# Create all needed folders
mkdir -p $DPATH $FPATH $LPATH

for ver in $VER
do
	# Create the folder needed
	rm -fr $LPATH/$ver/*
	mkdir -p $LPATH/$ver

	# Clean and copy all files needed
	cp -ft $LPATH/$ver $SPATH/submit$ver.sh $SPATH/$ver.sh $FILE

	# Move to the directory submit the code and return
	cd $LPATH/$ver
	bsub < submit$ver.sh
	cd ../../../
done