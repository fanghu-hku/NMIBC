#!/usr/bin/python
####################################################
#yanzeqin
#2022/4/12
#WGBS甲基化分析流程：QC、比对、去重、甲基化提取
#######################################################

#configfile: "config/WGBS.config.yaml"
SAMPLES = {}
with open(config["params"]["samples"], 'rb') as sinfo:
	for line in sinfo:
		parts = line.decode('utf-8').split()
		sample = parts[0]
		SAMPLES[sample] = [parts[1],parts[2]]

#print(SAMPLES)
rule all:
	input:expand(os.path.join(config["output"]["relative"],"{sample}", config["output"]["FilterCoverage"],"{sample}_bedGraph_Filtered.txt.gz"),sample=SAMPLES.keys())


##########################################01.QC-fastp##########################################################################
###############################################################################################################################

rule QC_fastp:
	input:
		raw_fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
		raw_fq2 = lambda wildcards: SAMPLES[wildcards.sample][1]
	output:
		fastp_fq1 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz"),
		fastp_fq2 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz"),
		fastp_html = os.path.join(config["output"]["relative"],"{sample}", config["output"]["fastp"], "{sample}.report.html"),
		fastp_json = os.path.join(config["output"]["relative"],"{sample}", config["output"]["fastp"], "{sample}.fastp.json")
	params:
		Sample = "{sample}",
		fastp = config["softwares"]["fastp"]
	shell:
		"""
		{params.fastp} -i {input.raw_fq1} -o {output.fastp_fq1} \
		-I {input.raw_fq2} -O {output.fastp_fq2} --compression 6 --report_title {params.Sample} \
		--json {output.fastp_json} --html {output.fastp_html} --detect_adapter_for_pe
		"""


rule bismark:
	input:
		fastp_fq1 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz"),
		fastp_fq2 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz")
	output:
		bam = os.path.join(config["output"]["relative"],"{sample}", config["output"]["bismark"], "{sample}.fastp_1_bismark_bt2_pe.bam")
	params:
		Sample = "{sample}",
		outdir = os.path.join(config["output"]["relative"],"{sample}", config["output"]["bismark"]),
		bismark = config["softwares"]["bismark"],
		bowtie2_path = config["softwares"]["bowtie2_path"],
		bismark_index = config["databases"]["bismark_index"]
	shell:
		"""
		{params.bismark} --nucleotide_coverage --path_to_bowtie2 {params.bowtie2_path} --parallel 4 --rg_id {params.Sample} --rg_sample {params.Sample} \
		--output_dir {params.outdir} {params.bismark_index} -1 {input.fastp_fq1} -2 {input.fastp_fq2}
		"""

rule deduplicate_bismark:
	input:
		bam = os.path.join(config["output"]["relative"],"{sample}", config["output"]["bismark"], "{sample}.fastp_1_bismark_bt2_pe.bam")
	output:
		dedup_bam = os.path.join(config["output"]["relative"],"{sample}", config["output"]["bismark"], "{sample}.deduplicated.bam")
	params:
		Sample = "{sample}",
		outdir = os.path.join(config["output"]["relative"],"{sample}", config["output"]["bismark"]),
		deduplicate_bismark = config["softwares"]["deduplicate_bismark"]
	shell:
		"""
		{params.deduplicate_bismark} --paired --bam --outfile {params.Sample} --output_dir {params.outdir} {input.bam}
		"""

rule bismark_methylation_extractor_1:
	input:
		dedup_bam = os.path.join(config["output"]["relative"],"{sample}", config["output"]["bismark"], "{sample}.deduplicated.bam")
	output:
		temp = os.path.join(config["output"]["relative"],"{sample}", config["output"]["methylation"],"bismark_methylation_extractor_1")
	params:
		Sample = "{sample}",
		metdir = os.path.join(config["output"]["relative"],"{sample}", config["output"]["methylation"]),
		extractor = config["softwares"]["bismark_methylation_extractor"]
	shell:
		"""
		{params.extractor} --paired-end --mbias_only --parallel 3 --output {params.metdir} {input.dedup_bam}
		touch {output.temp}
		"""

rule bismark_methylation_extractor_2:
	input:
		dedup_bam = os.path.join(config["output"]["relative"],"{sample}", config["output"]["bismark"], "{sample}.deduplicated.bam")
	output:
		temp = os.path.join(config["output"]["relative"],"{sample}", config["output"]["methylation"],"bismark_methylation_extractor_2")
	params:
		Sample = "{sample}",
		metdir = os.path.join(config["output"]["relative"],"{sample}", config["output"]["methylation"]),
		extractor = config["softwares"]["bismark_methylation_extractor"]
	shell:
		"""
		{params.extractor} --parallel 4 --ignore 5 --ignore_r2 5 --buffer_size 36G --output {params.metdir} \
		--paired-end --comprehensive --no_header --mbias_off --bedGraph --gzip {input.dedup_bam}
		touch {output.temp}
		"""

rule FilterCoverage:
	input:
		temp1 = os.path.join(config["output"]["relative"],"{sample}", config["output"]["methylation"],"bismark_methylation_extractor_1"),
		temp2 = os.path.join(config["output"]["relative"],"{sample}", config["output"]["methylation"],"bismark_methylation_extractor_2")
	output:
		bed = os.path.join(config["output"]["relative"],"{sample}", config["output"]["FilterCoverage"],"{sample}_bedGraph_Filtered.txt.gz")
	params:
		Sample = "{sample}",
		metdir = os.path.join(config["output"]["relative"],"{sample}", config["output"]["methylation"]),
		filt = os.path.join(config["output"]["relative"],"{sample}", config["output"]["FilterCoverage"]),
		python = config["softwares"]["python"],
		filter_script = config["scripts"]["filter_coverage"]
	shell:
		"""
		{params.python} {params.filter_script} --input {params.metdir} --output {params.filt} --basename {params.Sample} --coverage 5
		"""
