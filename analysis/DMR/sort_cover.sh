#!/bin/bash
# sort_cover.sh
# Process methylation data: remove chr prefix, adjust coordinates, calculate methylation ratio, sort
# Usage: sh sort_cover.sh <input_file> <output_file> <temp_file> <bedtools_sorted> <chr_filtered>

input_file=$1
output_file=$2
temp_file=$3
temp_f2=$4
temp_f3=$5

# bedtools sort -i ${input_file} > ${temp_f2}
cat ${temp_f2} | grep chr > ${temp_f3}

# Process data using awk
# Remove chr prefix from first column
# Subtract 1 from second column (0-based coordinates)
# Convert fourth column to methylation ratio
# Keep only first four columns

awk '{
    gsub(/^chr/, "", $1);
    $2 = $2 - 1;
    $4 = $4 / 100;
    print $1"\t"$2"\t"$3"\t"$4;
}' OFS="\t" ${temp_f3} > $temp_file

# Sort processed data and output to final file
sort -k1,1V -k2,2n $temp_file > $output_file

echo "finish 1"

# Clean up temporary files
rm $temp_file
