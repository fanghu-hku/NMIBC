#!/bin/bash
# sample.sh
# Classify Tumor/Normal type based on sample name
# Usage: sh sample.sh <input_file> <output_file>

input_file=$1
output_file=$2

# Read input file and process each line
cat $input_file | while IFS=$'\t' read -r name path
do
    if [[ $name == *"N"* ]]; then
        echo -e "Normal\t${name}\t${path}"
    elif [[ $name == *"T"* ]]; then
        echo -e "Tumor\t${name}\t${path}"
    else
        echo -e "None\t${name}\t${path}"
    fi
done > $output_file

# View output results
head $output_file
