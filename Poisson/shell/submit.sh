#!/bin/sh

# Meta vars
DATA=results/
LOGS=logs/
FIGPTH=figures/


# Define files needed by the execution in all tests
FILE="jacobiSolver.bin"

# Define which shell script will be executed
VER="mflop"

# Do stuff
mkdir -p $DATA $LOGS $FIGPTH
for ver in $VER
do
	# Copy all files
	mkdir $VER -p
	cp submit$ver.sh $ver/submit$ver.sh -f
	cp $ver.sh $ver/$ver.sh -f

	for f in $FILE
	do
		cp $f $ver/$f -f
	done

	# Move to the directory submit the code and return
	cd $ver/
	bsub < submit$ver.sh
	cd ../
done
