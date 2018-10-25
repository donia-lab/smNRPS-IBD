#!/bin/bash
# To output the total reads in a text file for each sample 

INPUT_FILE="$1"
OUTPUT_FILE="$2"


if [ ! -e $2 ]; then #if output doesn't exist
	result=$(grep -c '>' $INPUT_FILE) #find total reads
	echo "$result"> $2 #output number in an output file
fi
