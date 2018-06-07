#!/bin/sh

# Read in number of cores option
if [ -z "$1" ] ; then
	NUMCORES=4
else
	NUMCORES="$1"
fi

# Start interactive node
echo "Executing the command:"
echo "bsub -W 4:00 -q hpcint -app qrsh -env \"TERM\" -Is -n $NUMCORES \"reset; $SHELL -l\""
echo " "

bsub -J interactive_login -W 4:00 -q hpcint -app qrsh -env "TERM" -Is -n $NUMCORES "reset; $SHELL -l"

exit 0