#!/bin/sh
# Define meta variables

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
		./jacobiSolver.bin $met $n $n $n >> $met.dat
	done
done
