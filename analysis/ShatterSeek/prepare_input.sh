#!/bin/bash
# prepare_input.sh
# Prepare ShatterSeek input files
# Usage: sh prepare_input.sh <sample_name> <svaba_anno> <manta_anno> <jabba_seg> <output_dir>

SAMPLE=$1
SVABA_ANNO=$2
MANTA_ANNO=$3
JABBA_SEG=$4
OUTDIR=$5

mkdir -p ${OUTDIR}

# Generate CNV input file
CN_OUT=${OUTDIR}/cn.csv
echo -e "chrom\tstart\tend\tCN" > ${CN_OUT}
cat ${JABBA_SEG} | grep -v "track.name" | grep -v KI | grep -v GL | grep -vw Y | grep -wv M | \
    awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$6}' >> ${CN_OUT}

# Generate SV input file
SV_OUT=${OUTDIR}/sv.csv
echo -e "chrom1\tstart1\tchrom2\tend2\tSVtype\tstrand1\tstrand2" > ${SV_OUT}
cat ${SVABA_ANNO} ${MANTA_ANNO} | grep -v chrom1 | grep -vw "None" | grep -vw Y | grep -wv M | \
    awk -F "\t" '{print $1"\t"$2"\t"$8"\t"$9"\t"$16"\t"$6"\t"$13}' >> ${SV_OUT}

echo "Input files generated:"
echo "  CNV: ${CN_OUT}"
echo "  SV: ${SV_OUT}"
