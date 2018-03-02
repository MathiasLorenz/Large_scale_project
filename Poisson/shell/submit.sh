#!/bin/sh

# Run make command to ensure 
# Meta vars
DATA=results/
LOGS=logs/
FIGPTH=figures/
SHDIR=shell/


# Define files needed by the execution in all tests
FILE="jacobiSolver.bin"

# Define which shell script will be executed
VER="mflop"

# Do stuff
mkdir -p $DATA $LOGS $FIGPTH

mkdir $VER -p
for ver in $VER
do
	ls
	# Copy all files
	cp -ft $LOGS/$ver shell/submit$ver.sh shell/$ver.sh $FILE

	# Move to the directory submit the code and return
	cd $ver/
	bsub < submit$ver.sh
	cd ../
done
