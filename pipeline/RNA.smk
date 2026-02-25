#!/usr/bin/python
####################################################
#yanzeqin
#2022/4/12
#此流程步骤为lncRNA的基因融合检测
#######################################################

#configfile: "config.yaml"
SAMPLES = {}
with open(config["params"]["samples"], 'rb') as sinfo:
    for line in sinfo:
        parts = line.decode('utf-8').split()
        sample = parts[0]
        SAMPLES[sample] = [parts[1],parts[2]]

#print(SAMPLES)
rule all:
	input:expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["star_fusion"],"star-fusion.fusion_predictions.tsv"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["TE_salmon"],"quant.sf"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["trust4"],"{sample}_report.tsv"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken"],"{sample}_1.bracken"),sample=SAMPLES.keys()),
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["salmon"],"quant.sf"),sample=SAMPLES.keys())


##########################################01.QC-fastp##########################################################################
###############################################################################################################################
rule QC_fastp:
	input:
		raw_fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
		raw_fq2 = lambda wildcards: SAMPLES[wildcards.sample][1],
	output:
		fastp_fq1 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz"),
		fastp_fq2 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz"),
		fastp_html = os.path.join(config["output"]["relative"],"{sample}", config["output"]["fastp"], "{sample}.report.html"),
		fastp_json = os.path.join(config["output"]["relative"],"{sample}", config["output"]["fastp"], "{sample}.fastp.json"),
		fastp1_stat = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1_seqkit_stats.txt"),
		fastp2_stat = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2_seqkit_stats.txt"),
		raw1_stat = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.raw_1_seqkit_stats.txt"),
		raw2_stat = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.raw_2_seqkit_stats.txt")
	params:
		fastp = config["softwares"]["fastp"],
		seqkit = config["softwares"]["seqkit"],
		fastp_threads = config["fastp"]["threads"],
		f_adapter = config["fastp"]["f_adapter"],
		r_adapter = config["fastp"]["r_adapter"],
		quality = config["fastp"]["q"],
		un_qualified = config["fastp"]["u"],
		length_required = config["fastp"]["length_required"]
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "logs/{sample}.fastp.logs")
	benchmark:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "benchmarks/{sample}.fastp.benchmark.txt")
	shell:
		'''
		echo {input.raw_fq1}
		{params.fastp} -w {params.fastp_threads} \
		-i {input.raw_fq1}  -I {input.raw_fq2} \
		-q {params.quality} \
		-u 5 \
		-n 6 \
		-o {output.fastp_fq1} \
		-O {output.fastp_fq2} \
		-h {output.fastp_html} -j {output.fastp_json}
		{params.seqkit} stats -j 2 -a  {output.fastp_fq1} > {output.fastp1_stat}
		{params.seqkit} stats -j 2 -a  {output.fastp_fq2} > {output.fastp2_stat}
		{params.seqkit} stats -j 2 -a  {input.raw_fq1} > {output.raw1_stat}
		{params.seqkit} stats -j 2 -a  {input.raw_fq2} > {output.raw2_stat}
		'''





rule star_fusion:
	input:
		clean_fq1 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz"),
		clean_fq2 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz")
	output:
		star_bam = os.path.join(config["output"]["relative"],"{sample}",config["output"]["star_fusion"],"Aligned.out.bam"),
		star_abridged = os.path.join(config["output"]["relative"],"{sample}",config["output"]["star_fusion"],"star-fusion.fusion_predictions.abridged.tsv"),
		star_predictions=os.path.join(config["output"]["relative"],"{sample}",config["output"]["star_fusion"],"star-fusion.fusion_predictions.tsv")
	params:
		output_dir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["star_fusion"])
	log:
		os.path.join(config["output"]["relative"], "{sample}", config["output"]["star_fusion"], "logs/{sample}.star_fusion.logs")
	shell:
		'''
		conda init bash
		source /share/Data01/yanzeqin/software/snakemake/conda/etc/profile.d/conda.sh
		conda activate /share/Data01/yanzeqin/software/snakemake/conda/envs/star-fusion1.6.0
		STAR-Fusion --genome_lib_dir /datapool/yanzeqin/database/df_starfusion/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/  \
		--left_fq {input.clean_fq1}     \
		--right_fq {input.clean_fq2}    \
		--output_dir {params.output_dir}
		'''

rule TE_salmon:
	input:
		clean_fq1 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz"),
		clean_fq2 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz")
	output:
		output_quant = os.path.join(config["output"]["relative"],"{sample}",config["output"]["TE_salmon"],"quant.sf")
	params:
		output_dir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["TE_salmon"])
	shell:
		'''
		conda init bash
		source /share/Data01/yanzeqin/software/snakemake/conda/etc/profile.d/conda.sh
		conda activate /share/Data01/yanzeqin/software/snakemake/conda/envs/salmon
		salmon quant -i /datapool/zhuguanghui/multi.omics/yanzeqin/software/TE/REdiscoverTE/REdiscoverTE  \
		--thread 4 --seqBias --libType A  \
		-1 {input.clean_fq1}  \
		-2 {input.clean_fq2}   \
		-o {params.output_dir}
		'''

rule REdiscoverTE:
	input:
		expand(os.path.join(config["output"]["relative"],"{sample}",config["output"]["TE_salmon"],"quant.sf"),sample=SAMPLES.keys())
	output:
		sampleslist = os.path.join(config["output"]["relative"] +"salmon_result_20220818"),
		result = directory(config["output"]["relative"] +"/REdiscoverTE/")
	params:
		sam="sample"
	shell:
		"""
		echo -e "{params.sam}\tquant_sf_path" > {output.sampleslist}
		ls -1 {input}  |awk -F / '{{print $7"\t"$0}}' >> {output.sampleslist}
		conda init bash
		source /share/Data01/yanzeqin/software/snakemake/conda/etc/profile.d/conda.sh
		conda activate /share/Data01/yanzeqin/software/snakemake/conda/envs/salmon
		/share/Data01/yanzeqin/software/snakemake/conda/envs/r4/bin/Rscript /datapool/zhuguanghui/multi.omics/yanzeqin/software/TE/REdiscoverTE/rollup.R \
		--metadata={output.sampleslist} \
		--datadir=/datapool/zhuguanghui/multi.omics/yanzeqin/software/TE/REdiscoverTE/rollup_annotation/ \
		--nozero --threads=10 --assembly=hg38 \
		--outdir={output.result}
		"""


rule Salmon:
	input:
		clean_fq1 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz"),
		clean_fq2 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz")
	output:
		output_quant = os.path.join(config["output"]["relative"],"{sample}",config["output"]["salmon"],"quant.sf")
	params:
		output_dir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["salmon"])
	shell:
		"""
		/share/Data01/yanzeqin/software/snakemake/conda/envs/salmon/bin/salmon \
		quant -i /share/database/openData/GRCh38_GENCODE/transcripts_index_salmon \
		--gcBias -l A --thread 4 -1 {input.clean_fq1} -2 {input.clean_fq2} -o {params.output_dir}
		"""

rule trust4:
	input:
		clean_fq1 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz"),
		clean_fq2 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz")
	output:
		output_quant = os.path.join(config["output"]["relative"],"{sample}",config["output"]["trust4"],"{sample}_report.tsv")
	params:
		output_dir = os.path.join(config["output"]["relative"],"{sample}",config["output"]["trust4"]),
		name = os.path.join(config["output"]["relative"],"{sample}",config["output"]["trust4"],"{sample}")
	shell:
		'''
		mkdir -p {params.output_dir}
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/zhuguanghui/multi.omics/yzqsoftware/miniconda/envs/trust4
		/datapool/zhuguanghui/multi.omics/yzqsoftware/miniconda/envs/trust4/bin/run-trust4 -f /datapool/yanzeqin/database/trust4/TRUST4/human_IMGT+C.fa --ref /datapool/yanzeqin/database/trust4/TRUST4/human_IMGT+C.fa -1 {input.clean_fq1} -2 {input.clean_fq2} -o {params.name} -t 8 
		'''

rule karken:
	input:
		clean_fq1 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_1.fq.gz"),
		clean_fq2 = os.path.join(config["output"]["relative"], "{sample}", config["output"]["fastp"], "{sample}.fastp_2.fq.gz")
	output:
		kra=os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken"],"{sample}_1.kraken"),
		kreport=os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken"],"{sample}_1.kreport"),
		bra=os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken"],"{sample}_1.bracken")
	params:
		dir=os.path.join(config["output"]["relative"],"{sample}",config["output"]["kraken"])
	shell:
		'''
		mkdir -p {params.dir}
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/yanzeqin/software/anaconda/envs/kraken2
		kraken2 --paired --gzip-compressed --use-names --threads 8 --output {output.kra} --report {output.kreport} --db /datapool/zhuguanghui/multi.omics/stLFR/kraken2db   \
		{input.clean_fq1} {input.clean_fq2}
		bracken -d /datapool/zhuguanghui/multi.omics/stLFR/kraken2db -i {output.kreport} -o {output.bra} -t 8
		'''

