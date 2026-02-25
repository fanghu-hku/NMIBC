#!/bin/bash
# run_metilene_input.sh
# Generate metilene input file
# Usage: sh run_metilene_input.sh <tumor_name_file> <normal_name_file> <input_dir>

TUMOR_FILE=$1
NORMAL_FILE=$2
INPUT_DIR=$3

# Initialize variables (do not start with comma)
Input1=""
Input2=""

# Read tumor sample names and generate file paths
while IFS= read -r l1 || [ -n "$l1" ]; do
    file1="${INPUT_DIR}/${l1}/FilterCoverage/Tumor_${l1}_Coverage_Filtered.txt"
    if [ -f "$file1" ]; then
        Input1="${Input1},${file1}"
    else
        echo "Warning: $file1 does not exist."
    fi
done < $TUMOR_FILE

# Read normal sample names and generate file paths
while IFS= read -r l1 || [ -n "$l1" ]; do
    file2="${INPUT_DIR}/${l1}/FilterCoverage/Normal_${l1}_Coverage_Filtered.txt"
    if [ -f "$file2" ]; then
        Input2="${Input2},${file2}"
    else
        echo "Warning: $file2 does not exist."
    fi
done < $NORMAL_FILE

# Remove leading comma from strings
Input1="${Input1#,}"
Input2="${Input2#,}"

# Run metilene_input.pl to generate input file
metilene_input.pl \
	--in1 ${Input1} \
	--in2 ${Input2} \
	--h1 Tumor --h2 Normal

# Output: metilene_Tumor_Normal.input file
