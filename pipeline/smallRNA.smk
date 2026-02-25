#!/usr/bin/python
####################################################
# Author: yanzeqin
# Date: 2021-12-20
# small RNA pipeline: miRNA/piRNA QC, alignment, annotation, quantification
#######################################################


###########
# Libraries
###########
import pandas as pd

###############
# Configuration
###############
#configfile: "config/smallRNA.config.yaml"
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]


########################
# Samples and conditions
########################
SAMPLES = {}
with open(config["units"], 'rb') as sinfo:
	for line in sinfo:
		parts = line.decode('utf-8').split()
		sample = parts[0]
		SAMPLES[sample] = [parts[1],1]


###########################
# Input functions for rules
###########################
#################
# Desired outputs
#################


rule all:
	input:RESULT_DIR +"qc_report/"+"multiqc_report.html",RESULT_DIR +"qc_piRNAreport/"+"multiqc_report.html"


###############################
###        01.QC-fastp       ##
###############################


rule get_genome_annotation_gff:
	output:
		dna=config["databases"]["GRCh38"],
		trans=config["databases"]["transcripts_fa"],
		gentrome=config["databases"]["gentrome"],
		gtf=config["databases"]["gtf"],
		gff=config["databases"]["mirbase_gff"],
		bed=config["databases"]["mirbase_bed"],
		index2=multiext(config["databases"]["bowtie_index"],".1.ebwt",".2.ebwt",".3.ebwt",".4.ebwt",".rev.1.ebwt",".rev.2.ebwt"),
		index=directory(config["databases"]["salmon_index"])
	params:
		db_dir = config["databases"]["db_dir"],
		bowtie_build = config["softwares"]["bowtie_build"],
		salmon = config["softwares"]["salmon"]
	shell:
		"""
		wget -P {params.db_dir} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
		wget -P {params.db_dir} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz
		gzip -d GRCh38.primary_assembly.genome.fa.gz
		gzip -d gencode.v41.transcripts.fa.gz
		cat GRCh38.primary_assembly.genome.fa gencode.v41.transcripts.fa > gentrome.fasta
		{params.bowtie_build} GRCh38.primary_assembly.genome.fa bowtie_index
		{params.salmon} index -k 21 -t gencode.v41.transcripts.fa --gencode -i salmon_index
		wget -P {params.db_dir} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
		gzip -d gencode.v41.annotation.gtf.gz

		wget -P {params.db_dir} https://www.mirbase.org/ftp/22/genomes/hsa.gff3 --no-check-certificate
		cat hsa.gff3 | awk '{{if ($3=="miRNA") print}}' > mirbase.gff
		cat hsa.gff3 | awk 'BEGIN {{OFS = "\t"}};{{if ($3=="miRNA") print$1,$4-1,$5}}' > mirbase.bed
		#piRNA
		wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/bed/hsa.align.bed.gz
		wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/hsa.v3.0.fa.gz
		"""

rule fastp:
	input:
		lambda wildcards: SAMPLES[wildcards.sample][0]
	output:
		html = RESULT_DIR + "{sample}/" +"fastp/{sample}_fastp.html",
		zip = RESULT_DIR + "{sample}/" +"fastp/{sample}_fastp.fq.gz",
		json = RESULT_DIR + "{sample}/" +"fastp/{sample}_fastp.json"
	params:
		fastp = config["softwares"]["fastp"]
	log: RESULT_DIR + "{sample}/" + "fastp/{sample}.log"
	shell:
		"""
		{params.fastp} -w 4 \
		-i {input} \
		--adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
		--adapter_sequence_r2 GAACGACATGGCTACGATCCGACTT \
		-q 20 -u 5 --length_required 18 --length_limit 40 -n 2 \
		-o {output.zip} \
		-h {output.html} -j {output.json} 2> {log}
		"""

rule bowtie_align:
	input:
		reads= RESULT_DIR + "{sample}/" +"fastp/{sample}_fastp.fq.gz",
		index2=multiext(config["databases"]["bowtie_index"],".1.ebwt",".2.ebwt",".3.ebwt",".4.ebwt",".rev.1.ebwt",".rev.2.ebwt")
	output:
		temp(RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.sam")
	params:
		index=lambda w, input: os.path.commonprefix(input.index2).rstrip("."),
		bowtie = config["softwares"]["bowtie"]
	log: RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.log"
	shell:
		"""
		{params.bowtie} \
		--threads 8 \
		-n 1 \
		-l 20 \
		--best \
		-m 1 \
		--strata \
		--sam \
		{input.reads} \
		-x {params.index} \
		{output} 2>{log}
		"""

rule samtools_bam:
	input:
		RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.sam"
	output:
		temp(RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.bam")
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} view -bS {input} > {output}
		"""

rule samtools_sort:
	input:
		RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.bam"
	output:
		RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.sorted.bam"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} sort -o {output} {input}
		"""

rule samtools_index:
	input:
		RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.sorted.bam"
	output:
		RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.sorted.bam.bai"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} index {input}
		"""

rule samtools_flagstat:
	input:
		RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.sorted.bam"
	output:
		RESULT_DIR + "{sample}/" +"samtools/{sample}.flagstat"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} flagstat {input} > {output}
		"""

rule samtools_filter_mirbase:
	input:
		bam=RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.sorted.bam",
		bed=config["databases"]["mirbase_bed"]
	output:
		RESULT_DIR + "{sample}/" +"filtered_bam/{sample}.bam"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} view -b -h -L {input.bed} {input.bam}> {output}
		"""

rule samtools_filter_mirbase_index:
	input:
		RESULT_DIR + "{sample}/" +"filtered_bam/{sample}.bam"
	output:
		RESULT_DIR + "{sample}/" +"filtered_bam/{sample}.bam.bai"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} index {input}
		"""

rule bedtools_annotate:
	input:
		bam=RESULT_DIR + "{sample}/" +"filtered_bam/{sample}.bam",
		gff=config["databases"]["mirbase_gff"]
	output:
		RESULT_DIR + "{sample}/" +"annotation/{sample}.bed"
	shell:
		"""
		intersectBed -abam {input.bam} -b {input.gff} -r -wa -wb -bed > {output}
		"""

rule featureCounts:
	input:
		bam=RESULT_DIR + "{sample}/" +"filtered_bam/{sample}.bam",
		gff=config["databases"]["mirbase_gff"],
		fasta=config["databases"]["GRCh38"]
	output:
		multiext(RESULT_DIR + "{sample}/" +"featureCounts/{sample}",".featureCounts", ".featureCounts.summary")
	params:
		featureCounts = config["softwares"]["featureCounts"]
	log: RESULT_DIR + "{sample}/" +"featureCounts/{sample}.txt"
	shell:
		"""
		{params.featureCounts} \
		-T 2 \
		-t miRNA \
		-g Name \
		-G {input.fasta} \
		-F GFF \
		--fracOverlap 0.2 \
		-O \
		-s 0 \
		-M \
		-a \
		{input.gff} \
		-o {output[0]} {input.bam} 2> {log}
		"""

rule samtools_bam2fastq:
	input:
		RESULT_DIR + "{sample}/" + "filtered_bam/{sample}.bam"
	output:
		RESULT_DIR + "{sample}/" + "filtered_fastq/{sample}.fastq.gz"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} fastq {input} | gzip > {output}
		"""

rule salmon_quant:
	input:
		fastq=RESULT_DIR + "{sample}/" +"filtered_fastq/{sample}.fastq.gz",
	output:
		quant=RESULT_DIR + "{sample}/" +"salmon/{sample}_quant/quant.sf",
		result=directory(RESULT_DIR + "{sample}/" +"salmon/{sample}_quant")
	params:
		outdir=RESULT_DIR + "{sample}/" +"salmon/{sample}",
		index=config["databases"]["salmon_index"],
		salmon = config["softwares"]["salmon"]
	shell:
		"""
		{params.salmon} quant \
		-i {params.index} \
		-l A \
		-p 4 \
		-r {input.fastq} \
		--validateMappings \
		-o {params.outdir}_quant
		"""
rule multiqc:
	input:
		expand([RESULT_DIR + "{sample}/" +"fastp/{sample}_fastp.json",
			RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.log",
			RESULT_DIR + "{sample}/" +"samtools/{sample}.flagstat",
			RESULT_DIR + "{sample}/" +"featureCounts/{sample}.featureCounts.summary",
			RESULT_DIR + "{sample}/" +"salmon/{sample}_quant/"], sample=SAMPLES.keys())
	output:
		RESULT_DIR +"qc_report/"+"multiqc_report.html"
	params:
		out2 =  RESULT_DIR  +"qc_report/",
		multiqc = config["softwares"]["multiqc"]
	shell:
		"""
		{params.multiqc} {input} -n multiqc_report -f -q -o {params.out2}
		"""

rule samtools_filter_piRNA:
	input:
		bam=RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.sorted.bam",
		bed=config["databases"]["pirna_bed"]
	output:
		RESULT_DIR + "{sample}/" +"filtered_piRNAbam/{sample}.bam"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} view -b -h -L {input.bed} {input.bam}> {output}
		"""

rule samtools_filter_piRNA_index:
	input:
		RESULT_DIR + "{sample}/" +"filtered_piRNAbam/{sample}.bam"
	output:
		RESULT_DIR + "{sample}/" +"filtered_piRNAbam/{sample}.bam.bai"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} index {input}
		"""


rule featureCounts_piRNA:
	input:
		bam=RESULT_DIR + "{sample}/" +"filtered_piRNAbam/{sample}.bam",
		gff=config["databases"]["pirna_gff"],
		fasta=config["databases"]["GRCh38"]
	output:
		multiext(RESULT_DIR + "{sample}/" +"featureCounts_piRNA/{sample}",".featureCounts", ".featureCounts.summary")
	params:
		featureCounts = config["softwares"]["featureCounts"]
	log: RESULT_DIR + "{sample}/" +"featureCounts_piRNA/{sample}.txt"
	shell:
		"""
		{params.featureCounts} \
		-T 2 \
		-t piRNA \
		-g Name \
		-G {input.fasta} \
		-F GFF \
		--fracOverlap 0.2 \
		-O \
		-s 0 \
		-M \
		-a \
		{input.gff} \
		-o {output[0]} {input.bam} 2> {log}
		"""

rule samtools_bam2fastq_piRNA:
	input:
		RESULT_DIR + "{sample}/" + "filtered_piRNAbam/{sample}.bam"
	output:
		RESULT_DIR + "{sample}/" + "filtered_piRNAfastq/{sample}.fastq.gz"
	params:
		samtools = config["softwares"]["samtools"]
	shell:
		"""
		{params.samtools} fastq {input} | gzip > {output}
		"""

rule salmon_quant_piRNA:
	input:
		fastq=RESULT_DIR + "{sample}/" +"filtered_piRNAfastq/{sample}.fastq.gz",
	output:
		quant=RESULT_DIR + "{sample}/" +"piRNA_salmon/{sample}_quant/quant.sf",
		result=directory(RESULT_DIR + "{sample}/" +"piRNA_salmon/{sample}_quant")
	params:
		outdir=RESULT_DIR + "{sample}/" +"piRNA_salmon/{sample}",
		index=config["databases"]["salmon_index"],
		salmon = config["softwares"]["salmon"]
	shell:
		"""
		{params.salmon} quant \
		-i {params.index} \
		-l A \
		-p 4 \
		-r {input.fastq} \
		--validateMappings \
		-o {params.outdir}_quant
		"""

rule multiqc_piRNA:
	input:
		expand([RESULT_DIR + "{sample}/" +"fastp/{sample}_fastp.json",
			RESULT_DIR + "{sample}/" +"bowtie_align/{sample}.log",
			RESULT_DIR + "{sample}/" +"samtools/{sample}.flagstat",
			RESULT_DIR + "{sample}/" +"featureCounts_piRNA/{sample}.featureCounts.summary",
			RESULT_DIR + "{sample}/" +"piRNA_salmon/{sample}_quant/"], sample=SAMPLES.keys())
	output:
		RESULT_DIR +"qc_piRNAreport/"+"multiqc_report.html"
	params:
		out2 =  RESULT_DIR  +"qc_piRNAreport/",
		multiqc = config["softwares"]["multiqc"]
	shell:
		"""
		{params.multiqc} {input} -n multiqc_report -f -q -o {params.out2}
		"""
