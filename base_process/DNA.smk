####################################################
#yanzeqin
#2021/12/20
#此流程步骤为小RNA的质控，比对，注释，定量
#######################################################


###########
# Libraries
###########
import pandas as pd

###############
# Configuration
###############
#configfile: "/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Script/config/config.yaml" # where to find parameters

RESULT_DIR = config["output"]["relative"]
T="T"
N="B"

########################
# Samples and conditions
########################
SAMPLES = {}
with open(config["units"], 'rb') as sinfo:
        for line in sinfo:
                parts = line.decode('utf-8').split()
                sample = parts[0]
                SAMPLES[sample] = [parts[1],parts[2],parts[3],parts[4]]


###########################
# Input functions for rules
###########################
#################
# Desired outputs
#################


#print(SAMPLES)
rule all:
	input:
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["mutect2withpon"] + "/{sample}_Unfiltered.vcf"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["GetPileupSummaries"] + "/{sample}_"+N+"_PileupSummary.table"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["FilterMutect2WithPON"] + "/{sample}_Filtered.vcf"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["manta"] + "/results/variants/somaticSV.vcf.gz"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_BQSR.bam_CNVs"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+".seg"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["sequenza"] + "/result/{sample}_confints_CP.txt"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["GermlineSNPs"] + "/{sample}_"+N+".vcf.gz"],sample=SAMPLES.keys()),
		#RESULT_DIR +"/"+config["output"]["GermlineSNPs_db2"] + "/GermlineSNPs_filted.vcf",
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.bed"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["AA_cnvkit"] + "/{sample}_AA_CNV_SEEDS.bed"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/" + config["output"]["cnvkit_GISTIC"] + "/all_lesions.conf_90.txt"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["AmpliconArchitect2"] + "/{sample}_summary.txt"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["kraken"] + "/{sample}_"+T+".bracken"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["svaba"] + "/{sample}.svaba.somatic.sv.vcf"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["strelka"] + "/results/variants/somatic.snvs.vcf.gz"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["AnnotSV"] + "/{sample}_manta.tsv"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["AnnotSV"] + "/{sample}_svaba.tsv"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_clean.denoisedCR.tsv"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"] + "/{sample}_SURVIVOR.anno.txt"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"] + "/{sample}_SURVIVOR.tsv"],sample=SAMPLES.keys()),
		expand([RESULT_DIR + "/Result/sample_bamdepth/{sample}_"+T+"_bam.depth.gz"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["smoove"] + "/{sample}-smoove.genotyped.vcf.gz"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/{sample}/" + config["output"]["JaBbA"] + "/jabba.simple.cnv.vcf"],sample=SAMPLES.keys()),
		#expand([RESULT_DIR + "/" + config["output"]["GISTIC"] + "/all_lesions.conf_90.txt"],sample=SAMPLES.keys())
############################
#####01.QC-fastp############
############################
rule fastp_T:
	input:
		rawfq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
		rawfq2 = lambda wildcards: SAMPLES[wildcards.sample][1]
	output:
		read1 =  temp(RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"_fastp1.fq.gz"),
		read2 =  temp(RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"_fastp2.fq.gz"),
		fastp_html = RESULT_DIR + "/{sample}/" + config["output"]["fastp"] +  "/{sample}_"+T+".report.html",
		fastp_json = RESULT_DIR + "/{sample}/" + config["output"]["fastp"] +  "/{sample}."+T+".fastp.json"
	params:
		v1 = "{sample}_"+T,
	log: RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"fastp.log"
	shell:
                r'''
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/fastp \
		-i {input.rawfq1} -o {output.read1} -I {input.rawfq2} -O {output.read2} \
		--json {output.fastp_json} --html {output.fastp_html} \
		--report_title {params.v1} \
		--detect_adapter_for_pe --compression 6 --cut_front --cut_tail	 2> {log}
                '''

rule BWA_T:
	input:
		read1= RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"_fastp1.fq.gz",
		read2= RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"_fastp2.fq.gz"
	output:
		bam1= temp(RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}_"+T+"_Sorted.bam")
	params:
		rg=r"@RG\tID:{sample}_"+T+r"\tSM:{sample}_"+T+r"\tLB:{sample}_"+T+r"\tPU:{sample}_"+T+r"\tPL:ILLUMINA",
		GRCh38= config["databases"]["GRCh38"],
		v1 = "{sample}_"+T
	log: RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}_"+T+".log"
	shell:
		r'''
		/share/Data01/yanzeqin/software/snakemake/conda/bin/bwa \
		mem -M -t 8 -R "{params.rg}" {params.GRCh38} {input.read1} {input.read2}|/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/samtools \
		sort -m 4G -l 6 -O BAM -T {params.v1} --threads 6 -o {output.bam1} 2> {log}
		'''



rule MarkDuplicates_T:
	input:
		bam1= RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}_"+T+"_Sorted.bam"
	output:
		bam2= temp(RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_MD.bam"),
		bam3= temp(RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_MD_Sorted.bam"),
		metrics= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_dup_metrics.txt"
	params:
		v1 = "{sample}_"+T,
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"+T
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_MarkDuplicates.log"
	shell:
		'''
		mkdir -p {params.TempDir}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx100G -Xms16G" MarkDuplicates -I {input.bam1} -O {output.bam2} -M {output.metrics} --TMP_DIR {params.TempDir} 2> {log}
		echo -e "Gatk MarkDuplicates done!" >> {log}
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/samtools sort -l 6 -O bam -T {params.v1} --threads 8 -o {output.bam3} {output.bam2} 2>> {log}
		rm -rf {params.TempDir}
		'''

rule BQSR_T:
	input:bam= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_MD_Sorted.bam"
	output:table= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+".table" 
	params:
		GRCh38= config["databases"]["GRCh38"],
		dbSNP= config["databases"]["dbSNP"],
		G1000= config["databases"]["G1000"],
		ClinVar= config["databases"]["ClinVar"],
		COSMIC1= config["databases"]["COSMIC1"],
		COSMIC2= config["databases"]["COSMIC2"],
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"+T
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSRtable.log"
	shell:
		'''
		mkdir -p {params.TempDir}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx80G -Xms16G" BaseRecalibrator -I {input.bam} -R {params.GRCh38} \
		--known-sites {params.dbSNP} --known-sites {params.G1000} --known-sites {params.ClinVar} --known-sites {params.COSMIC1} --known-sites {params.COSMIC2} \
		-O {output.table} --bqsr-baq-gap-open-penalty 30 --tmp-dir {params.TempDir} 2> {log}
		rm -rf {params.TempDir}
		'''

rule ApplyBQSR_T:
	input:
		bam= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_MD_Sorted.bam",
		table= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+".table"	
	output:
		bam1= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	params:
		GRCh38= config["databases"]["GRCh38"],
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"+T
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_ApplyBQSR.log"
	shell:
		'''
		mkdir -p {params.TempDir}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx80G -Xms16G" ApplyBQSR -R {params.GRCh38} -I {input.bam} --bqsr-recal-file {input.table} -O {output.bam1} --tmp-dir {params.TempDir} 2> {log}
		rm -rf {params.TempDir}
		'''


rule fastp_N:
	input:
		rawfq1 = lambda wildcards: SAMPLES[wildcards.sample][2],
		rawfq2 = lambda wildcards: SAMPLES[wildcards.sample][3]
	output:
		read1 =  temp(RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+N+"_fastp1.fq.gz"),
		read2 =  temp(RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+N+"_fastp2.fq.gz"),
		fastp_html = RESULT_DIR + "/{sample}/" + config["output"]["fastp"] +  "/{sample}_"+N+".report.html",
		fastp_json = RESULT_DIR + "/{sample}/" + config["output"]["fastp"] +  "/{sample}."+N+".fastp.json"
	params:
		v1 = "{sample}_"+N,
	log: RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+N+"fastp.log"
	shell:
		'''
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/fastp \
		-i {input.rawfq1} -o {output.read1} -I {input.rawfq2} -O {output.read2} \
		--json {output.fastp_json} --html {output.fastp_html} \
		--report_title {params.v1} \
		--detect_adapter_for_pe --compression 6 --cut_front --cut_tail   2> {log}
		'''

rule BWA_N:
	input:
		read1= RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+N+"_fastp1.fq.gz",
		read2= RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+N+"_fastp2.fq.gz"
	output:
		bam1= temp(RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}_"+N+"_Sorted.bam")
	params:
		rg=r"@RG\tID:{sample}_"+N+r"\tSM:{sample}_"+N+r"\tLB:{sample}_"+N+r"\tPU:{sample}_"+N+r"\tPL:ILLUMINA",
		GRCh38= config["databases"]["GRCh38"],
		v1 = "{sample}_"+N
	log: RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}_"+N+".log"
	shell:
		r'''
		/share/Data01/yanzeqin/software/snakemake/conda/bin/bwa \
		mem -M -t 8 -R "{params.rg}" {params.GRCh38} {input.read1} {input.read2} |/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/samtools \
		sort -m 4G -l 6 -O BAM -T {params.v1} --threads 6 -o {output.bam1} 2> {log}
                '''



rule MarkDuplicates_N:
	input:
		bam1= RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}_"+N+"_Sorted.bam"
	output:
		bam2= temp(RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_MD.bam"),
		bam3= temp(RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_MD_Sorted.bam"),
		metrics= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_dup_metrics.txt"
	params:
		v1 = "{sample}_"+N,
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"+N
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_MarkDuplicates.log"
	shell:
		'''
		mkdir -p {params.TempDir}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx100G -Xms16G" MarkDuplicates -I {input.bam1} -O {output.bam2} -M {output.metrics} --TMP_DIR {params.TempDir} 2> {log}
		echo -e "Gatk MarkDuplicates done!" >> {log}
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/samtools sort -l 6 -O bam -T {params.v1} --threads 8 -o {output.bam3} {output.bam2} 2>> {log}
		rm -rf {params.TempDir}
		'''


rule BQSR_N:
	input:bam= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_MD_Sorted.bam"
	output:table= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+".table"
	params:
		GRCh38= config["databases"]["GRCh38"],
		dbSNP= config["databases"]["dbSNP"],
		G1000= config["databases"]["G1000"],
		ClinVar= config["databases"]["ClinVar"],
		COSMIC1= config["databases"]["COSMIC1"],
		COSMIC2= config["databases"]["COSMIC2"],
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"+N
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSRtable.log"
	shell:
		'''
                mkdir -p {params.TempDir}
                /share/Data01/pengguoyu/bin/gatk --java-options "-Xmx80G -Xms16G" BaseRecalibrator -I {input.bam} -R {params.GRCh38} \
                --known-sites {params.dbSNP} --known-sites {params.G1000} --known-sites {params.ClinVar} --known-sites {params.COSMIC1} --known-sites {params.COSMIC2} \
                -O {output.table} --bqsr-baq-gap-open-penalty 30 --tmp-dir {params.TempDir} 2> {log}
                rm -rf {params.TempDir}
		'''

rule ApplyBQSR_N:
	input:
		bam= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_MD_Sorted.bam",
		table= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+".table"
	output:
		bam1= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam"
	params:
		GRCh38= config["databases"]["GRCh38"],
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"+N
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_ApplyBQSR.log"
	shell:
		'''
		mkdir -p {params.TempDir}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx80G -Xms16G" ApplyBQSR -R {params.GRCh38} -I {input.bam} --bqsr-recal-file {input.table} -O {output.bam1} --tmp-dir {params.TempDir}  2> {log}
		rm -rf {params.TempDir}
		'''

rule samtools_depth:
	input:
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		depthN=RESULT_DIR + "/Result/sample_bamdepth/{sample}_"+N+"_bam.depth.gz",
		depthT=RESULT_DIR + "/Result/sample_bamdepth/{sample}_"+T+"_bam.depth.gz"
	
	shell:
		'''
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/samtools depth -aa {input.normal} > {output.depthN}
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/samtools depth -aa {input.tumor} > {output.depthT}
		'''

##########################################################################################################################################
rule Mutect2WithPON:
	input:
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		vcf= RESULT_DIR + "/{sample}/" + config["output"]["mutect2withpon"] + "/{sample}_Unfiltered.vcf", 
		fif2= RESULT_DIR + "/{sample}/" + config["output"]["mutect2withpon"] + "/{sample}_F1R2.tar.gz"
	params:
		pon=config["databases"]["pon"],
		GRCh38= config["databases"]["GRCh38"],
                GnomdAD= config["databases"]["GnomdAD_AF"],
		temp= RESULT_DIR + "/{sample}/" + config["output"]["mutect2withpon"] + "/TEMP",
		normal= "{sample}_"+N
	shell:
		'''
		mkdir -p {params.temp}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" Mutect2 --tmp-dir {params.temp} -R {params.GRCh38} -I {input.tumor} -I {input.normal} --normal-sample {params.normal} --germline-resource {params.GnomdAD} --panel-of-normals {params.pon} --f1r2-tar-gz {output.fif2} -O {output.vcf}
		rm -rf {params.temp}
		'''



rule GetPileupSummaries:
	input:
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		Ttable= RESULT_DIR + "/{sample}/" + config["output"]["GetPileupSummaries"] + "/{sample}_"+T+"_PileupSummary.table",
		Ntable= RESULT_DIR + "/{sample}/" + config["output"]["GetPileupSummaries"] + "/{sample}_"+N+"_PileupSummary.table"
		
	params:
		Case="{sample}",
                Normal="{sample}_"+N,
                Tumor="{sample}_"+T,
		TmpDir=RESULT_DIR + "/{sample}/TEMP",
		outDir=RESULT_DIR + "/{sample}/" + config["output"]["GetPileupSummaries"],
		GRCh38="/share/database/openData/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa",
		GRCh38Dict="/share/database/openData/GRCh38_GENCODE/GRCh38.primary_assembly.genome.dict",
		Gnomad="/share/database/openData/GATK_BestPracticesResource/somatic-hg38/af-only-gnomad.hg38.subset.vcf.gz",
		IntervalDir="/share/database/openData/GRCh38_GENCODE/splited_intervals",
		Intervals="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFPE2/Intervals"
	shell:
		'''
		mkdir -p {params.TmpDir}
		sh /datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFPE2/gem.sh {input.tumor} {input.normal} {params.outDir} {output.Ttable} {output.Ntable} {params.TmpDir} {params.Case} {params.Normal} {params.Tumor}
		rm -rf {params.TmpDir}
		'''


rule CalculateContamination:
	input:
		Ttable= RESULT_DIR + "/{sample}/" + config["output"]["GetPileupSummaries"] + "/{sample}_"+T+"_PileupSummary.table",
                Ntable= RESULT_DIR + "/{sample}/" + config["output"]["GetPileupSummaries"] + "/{sample}_"+N+"_PileupSummary.table",
	output:
		con= RESULT_DIR + "/{sample}/" + config["output"]["FilterMutect2WithPON"] + "/{sample}_Contamination.table"
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" CalculateContamination -I {input.Ttable} --matched-normal {input.Ntable} -O {output.con}
		'''


rule FilterMutect2WithPON:
	input:
		fif2= RESULT_DIR + "/{sample}/" + config["output"]["mutect2withpon"] + "/{sample}_F1R2.tar.gz",
		vcf= RESULT_DIR + "/{sample}/" + config["output"]["mutect2withpon"] + "/{sample}_Unfiltered.vcf",
		con= RESULT_DIR + "/{sample}/" + config["output"]["FilterMutect2WithPON"] + "/{sample}_Contamination.table",
		Ttable= RESULT_DIR + "/{sample}/" + config["output"]["GetPileupSummaries"] + "/{sample}_"+T+"_PileupSummary.table",
		Ntable= RESULT_DIR + "/{sample}/" + config["output"]["GetPileupSummaries"] + "/{sample}_"+N+"_PileupSummary.table"
	output:
		fileVcf= RESULT_DIR + "/{sample}/" + config["output"]["FilterMutect2WithPON"] + "/{sample}.vcf.gz",
		ReadModel= RESULT_DIR + "/{sample}/" + config["output"]["FilterMutect2WithPON"] + "/{sample}_ReadOrientModel.tar.gz",
		filted= RESULT_DIR + "/{sample}/" + config["output"]["FilterMutect2WithPON"] + "/{sample}_Filtered.vcf"
	params:
		fileVcf= RESULT_DIR + "/{sample}/" + config["output"]["FilterMutect2WithPON"] + "/{sample}.vcf",
		stats=  RESULT_DIR + "/{sample}/" + config["output"]["mutect2withpon"] + "/{sample}_Unfiltered.vcf.stats",
		Case="{sample}",
		Normal="{sample}_"+N,
		Tumor="{sample}_"+T
	shell:
		'''
		sh /datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFPE2/file.sh {input.fif2} {input.Ttable} {input.Ntable} {input.vcf} {params.stats} {output.ReadModel} {input.con} {params.fileVcf} {params.Case} {params.Normal} {params.Tumor} {output.fileVcf} {output.filted}
		'''

rule freec_config:
	input:
		con= RESULT_DIR + "/{sample}/" + config["output"]["FilterMutect2WithPON"] + "/{sample}_Contamination.table",
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		freectxt=RESULT_DIR + "/{sample}/{sample}_freec.txt"
	params:
		clinical="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Script/clinical_raw.csv",
		Case="{sample}",
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"]
	shell:
		'''
		sh /datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Script/get_config.sh {params.Case} {input.tumor} {input.normal} {params.outputDir} {params.clinical} {input.con} {output.freectxt}
		'''

rule freec:
	input:
		freectxt=RESULT_DIR + "/{sample}/{sample}_freec.txt"
	output:
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_BQSR.bam_CNVs",
		cpn= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_BQSR.bam_sample.cpn",
		ratio= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_BQSR.bam_ratio.txt"
	params:
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"]
	shell:
		'''
		mkdir -p {params.outputDir}
		/share/Data01/pengguoyu/bin/freec -conf {input.freectxt}
		'''


rule cnv_bed:
	input:
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_BQSR.bam_CNVs",
		ratio= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_BQSR.bam_ratio.txt"
	output:
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_CNVs.bed",
		ratio= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_ratio.bed"
	shell:
		'''
		python /datapool/zhuguanghui/multi.omics/stLFR/py_script/cnv_bed.py {input.outputDir} {output.outputDir}
		echo 1
		perl /share/Data01/pengguoyu/App/FREEC-11.6/scripts/freec2bed.pl -f {input.ratio} -p 2 > {output.ratio}
		'''

rule GISTIC_Seg:
	input:
		cpn= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_BQSR.bam_sample.cpn",
		ratio= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+"_ratio.bed"
	output:
		seg= RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+".seg"
	params:
		tumor="{sample}"+T
	shell:
		'''
		python /datapool/zhuguanghui/multi.omics/stLFR/py_script/GISTIC_Seg.py --bed {input.ratio} --cpn {input.cpn} --seg {output.seg} --sample {params.tumor}
		'''

rule GISTIC:
	input:expand([RESULT_DIR + "/{sample}/" + config["output"]["FREEC"] + "/{sample}_"+T+".seg"], sample=SAMPLES.keys())
	output:
		primaryseg= RESULT_DIR + "/" +config["output"]["GISTIC"] + "/primary_tumor.seg",
		result= RESULT_DIR + "/" +config["output"]["GISTIC"] + "/all_lesions.conf_90.txt"
	params:
		gisticdir= RESULT_DIR + "/" +config["output"]["GISTIC"]
	shell:
		'''
		mkdir -p {params.gisticdir}
		cat {input} > {output.primaryseg}
		/share/Data01/pengguoyu/bin/gistic2 -b {params.gisticdir} -seg {output.primaryseg} -refgene /share/Data01/pengguoyu/App/gistic2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -broad 1 -brlen 0.5 -conf 0.95 -maxseg 10000 -qvt 0.1 -ta 0.21 -td 0.25 -genegistic 1
		'''


rule manta:
	input:
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		conflow= RESULT_DIR + "/{sample}/" + config["output"]["manta"] + "/runWorkflow.py",
		result= protected(RESULT_DIR + "/{sample}/" + config["output"]["manta"] + "/results/variants/somaticSV.vcf.gz"),
		re= RESULT_DIR + "/{sample}/" + config["output"]["manta"] + "/results/variants/candidateSmallIndels.vcf.gz"
	params:
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["manta"]
	shell:
		'''
		rm -rf {params.outputDir}
		mkdir -p {params.outputDir}
		/share/Data01/pengguoyu/App/manta-1.6.0.centos6_x86_64/bin/configManta.py --normalBam {input.normal} --tumorBam {input.tumor} \
		--referenceFasta /datapool/zhuguanghui/multi.omics/stLFR/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa \
		--callRegions /datapool/zhuguanghui/multi.omics/stLFR/GRCh38_GENCODE/GRCh38.primary_assembly.genome.CommonChr.bed.gz \
		--runDir {params.outputDir}
		{output.conflow}
		'''

rule sequenza:
	input:
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		sample_seqz= RESULT_DIR + "/{sample}/" + config["output"]["sequenza"] + "/{sample}_sample.seqz.gz"
	params:
		hg38="/share/database/openData/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa",
		hg_wig="/share/service09/yanzeqin/yanzeqin/stLFR/sequenza/hg38.Base.wig.gz",
		id="{sample}",
		outfile= RESULT_DIR + "/{sample}/" + config["output"]["sequenza"] + "/result"
	shell:
		'''
		conda init bash
                source /share/Data01/yanzeqin/software/snakemake/conda/etc/profile.d/conda.sh
                conda activate /share/Data01/yanzeqin/software/snakemake/conda/envs/sequenza/
		/share/Data01/yanzeqin/software/snakemake/conda/envs/sequenza/bin/sequenza-utils bam2seqz -n {input.normal} -t {input.tumor} --fasta {params.hg38} -gc {params.hg_wig} -o {output.sample_seqz}
		'''

rule sequenza2:
	input:
		sample_seqz= RESULT_DIR + "/{sample}/" + config["output"]["sequenza"] + "/{sample}_sample.seqz.gz"
	output:
		samll= RESULT_DIR + "/{sample}/" + config["output"]["sequenza"] + "/{sample}_small.seqz.gz",
		result= RESULT_DIR + "/{sample}/" + config["output"]["sequenza"] + "/result/{sample}_confints_CP.txt"
	params:
		id="{sample}",
		outfile= RESULT_DIR + "/{sample}/" + config["output"]["sequenza"] + "/result"
	shell:
		'''
		mkdir -p {params.outfile}
		conda init bash
                source /share/Data01/yanzeqin/software/snakemake/conda/etc/profile.d/conda.sh
                conda activate /share/Data01/yanzeqin/software/snakemake/conda/envs/sequenza
		/share/Data01/yanzeqin/software/snakemake/conda/envs/sequenza/bin/sequenza-utils seqz_binning --seqz {input.sample_seqz} -w 50 -o {output.samll}
		/share/Data01/yanzeqin/software/snakemake/conda/envs/sequenza/bin/Rscript /datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Rscript/step2_sequenza.R {params.id} {output.samll} {params.outfile}
		'''

rule svaba:
	input:
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		result= RESULT_DIR + "/{sample}/" + config["output"]["svaba"] + "/{sample}.svaba.somatic.sv.vcf"
	params:
		outfile= RESULT_DIR + "/{sample}/" + config["output"]["svaba"],
		id="{sample}"
	shell:
		'''
		mkdir -p {params.outfile}
		cd {params.outfile}
		/datapool/yanzeqin/software/miniforge3/envs/svaba/bin/svaba run -t {input.tumor} -n {input.normal} \
		-G /datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa -a {params.id} -p 4
		'''

rule strelka:
	input:
		re= RESULT_DIR + "/{sample}/" + config["output"]["manta"] + "/results/variants/candidateSmallIndels.vcf.gz",
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		re = RESULT_DIR + "/{sample}/" + config["output"]["strelka"] + "/runWorkflow.py",
		result= RESULT_DIR + "/{sample}/" + config["output"]["strelka"] + "/results/variants/somatic.snvs.vcf.gz"
	params:
		outfile= RESULT_DIR + "/{sample}/" + config["output"]["strelka"]
	shell:
		'''
		mkdir -p {params.outfile}
		cd {params.outfile}
		conda init bash
                source /share/Data01/yanzeqin/software/snakemake/conda/etc/profile.d/conda.sh
                conda activate /share/Data01/yanzeqin/software/snakemake/conda/envs/sequenza
		/share/Data01/yanzeqin/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam {input.normal} --tumorBam {input.tumor} --referenceFasta /datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa --runDir {params.outfile}
		/usr/bin/python2 {output.re} -m local -j 4
		'''

rule AnnotSV_manta:
	input:
		result= RESULT_DIR + "/{sample}/" + config["output"]["manta"] + "/results/variants/somaticSV.vcf.gz"
	output:
	        re= RESULT_DIR + "/{sample}/" + config["output"]["AnnotSV"] + "/{sample}_manta.tsv"
	log: 
		RESULT_DIR + "/{sample}/" + config["output"]["AnnotSV"] + "/log"
	params:
		outname= "{sample}_manta",
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["AnnotSV"] 
	shell:
		'''
		conda init bash
		source /datapool/yanzeqin/software/miniforge3/etc/profile.d/conda.sh
		conda activate /datapool/yanzeqin/software/miniforge3/envs/AnnotSV
		echo "1"
		/datapool/yanzeqin/software/AnnotSV/AnnotSV_2.2/bin/AnnotSV/AnnotSV.tcl -SVinputFile {input.result} -outputDir {params.outputDir} -outputFile {params.outname} -genomeBuild GRCh38
		'''

rule AnnotSV_svaba:
	input:
		result= RESULT_DIR + "/{sample}/" + config["output"]["svaba"] + "/{sample}.svaba.somatic.sv.vcf"
	output:
		re= RESULT_DIR + "/{sample}/" + config["output"]["AnnotSV"] + "/{sample}_svaba.tsv"
	log:
		RESULT_DIR + "/{sample}/" + config["output"]["AnnotSV"] + "/svaba.log"
	params:
		outname= "{sample}_svaba",
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["AnnotSV"]
	shell:
		'''
		conda init bash
		source /datapool/yanzeqin/software/miniforge3/etc/profile.d/conda.sh
		conda activate /datapool/yanzeqin/software/miniforge3/envs/AnnotSV
		echo "1"
		/datapool/yanzeqin/software/AnnotSV/AnnotSV_2.2/bin/AnnotSV/AnnotSV.tcl -SVinputFile {input.result} -outputDir {params.outputDir} -outputFile {params.outname} -genomeBuild GRCh38
		'''

rule neoSV:
	input:
		re_manta= RESULT_DIR + "/{sample}/" + config["output"]["manta"] + "/results/variants/somaticSV.vcf.gz",
		re_svaba= RESULT_DIR + "/{sample}/" + config["output"]["svaba"] + "/{sample}.svaba.somatic.sv.vcf"
		
	output:
		temp = RESULT_DIR + "/{sample}/" + config["output"]["neoSV"] + "/{sample}_manta.vcf",
		re1 = RESULT_DIR + "/{sample}/" + config["output"]["neoSV"] + "/{sample}_manta.anno.txt",
		re2 = RESULT_DIR + "/{sample}/" + config["output"]["neoSV"] + "/{sample}_svaba.anno.txt",
	params:
		name1 = "{sample}_manta",
		name2 = "{sample}_svaba",
		outputdir = RESULT_DIR + "/{sample}/" + config["output"]["neoSV"]
	priority: 50
	shell:
		'''
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/zhuguanghui/multi.omics/yzqsoftware/miniconda/envs/neosv
		zcat {input.re_manta} > {output.temp}
		neosv -sf {output.temp} -o {params.outputdir} -p {params.name1} -r 95  --anno-only
		neosv -sf {input.re_svaba} -o {params.outputdir} -p {params.name2} -r 95  --anno-only
		
		'''


rule neo_mergeSV:
	input:
		re_del= RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"] + "/{sample}_SURVIVOR.vcf"
	output:
		re1 = RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"] + "/{sample}_SURVIVOR.anno.txt"
	params:
		name1 = "{sample}_SURVIVOR",
		outputdir = RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"]
	priority: 50
	shell:
		'''
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/zhuguanghui/multi.omics/yzqsoftware/miniconda/envs/neosv
		neosv -sf {input.re_del} -o {params.outputdir} -p {params.name1} -r 95  --anno-only
		'''


rule annotSV_mergeSV:
	input:
		re_del= RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"] + "/{sample}_SURVIVOR.vcf"

	output:
		re1= RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"] + "/{sample}_SURVIVOR.tsv"
	params:
		temp1= RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"] + "/temp1.vcf",
		outname1= "{sample}_SURVIVOR",
		outname2= "{sample}_otherAnnotSV",
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["mergeSV"]
		
	priority: 20
	shell:
		'''
		conda init bash
		source /datapool/yanzeqin/software/miniforge3/etc/profile.d/conda.sh
		conda activate /datapool/yanzeqin/software/miniforge3/envs/AnnotSV
		echo "1"
		/datapool/yanzeqin/software/AnnotSV/AnnotSV_2.2/bin/AnnotSV/AnnotSV.tcl -SVinputFile {input.re_del} -outputDir {params.outputDir} -outputFile {params.outname1} -genomeBuild GRCh38
		'''
	

rule GATK_GermlineSNPs2:
	input:
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam"
	output:
		vcf=RESULT_DIR + "/{sample}/" + config["output"]["GermlineSNPs"] + "/{sample}_"+N+".vcf.gz"
	params:
		GRCh38= config["databases"]["GRCh38"]
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" HaplotypeCaller -R {params.GRCh38} \
		-I {input.normal} \
		-O {output.vcf} \
		-contamination 0 -ERC GVCF
		'''

rule germlineSNPs_DB1:
	input:
		expand([RESULT_DIR + "/{sample}/" + config["output"]["GermlineSNPs"] + "/{sample}_"+N+".vcf.gz"],sample=SAMPLES.keys())
	output:
		dir=directory(RESULT_DIR +"/"+ config["output"]["GermlineSNPs_db"])
	params:
		V= " -V ".join(expand([RESULT_DIR + "/{sample}/" + config["output"]["GermlineSNPs"] + "/{sample}_"+N+".vcf.gz"],sample=SAMPLES.keys())),
		GRCh38= config["databases"]["GRCh38"],
		temp= RESULT_DIR  + "/TEMP_PON",
		intervals= config["databases"]["intervals"]
	shell:
		'''
		mkdir -p {params.temp}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx92G -Xms16G" GenomicsDBImport  --tmp-dir {params.temp} --genomicsdb-workspace-path {output.dir} -L {params.intervals} -R {params.GRCh38} -V {params.V} --batch-size 6 --genomicsdb-shared-posixfs-optimizations true
		rm -rf {params.temp}
		'''

rule GenotypeGVCFs:
	input:  ancient(RESULT_DIR +"/" +config["output"]["GermlineSNPs_db"])
	output:
		vcf=RESULT_DIR +"/"+config["output"]["GermlineSNPs_db2"] + "/GermlineSNPs_db.vcf"
	params:
		GRCh38= config["databases"]["GRCh38"],
		GnomdAD= config["databases"]["GnomdAD_AF"],
		temp= RESULT_DIR + "/" + config["output"]["GermlineSNPs_db2"]
	shell:
		'''
		mkdir -p {params.temp}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" GenotypeGVCFs \
		-R {params.GRCh38} \
		-V gendb://{input} \
		-O {output.vcf}
		'''

rule VariantFiltration:
	input:
		vcf=RESULT_DIR + "/"+ config["output"]["GermlineSNPs_db2"] + "/GermlineSNPs_db.vcf"
	output:
		vcf=RESULT_DIR +"/"+ config["output"]["GermlineSNPs_db2"] + "/GermlineSNPs_filted.vcf"
	params:
		GRCh38= config["databases"]["GRCh38"]
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx92G -Xms16G" VariantFiltration \
		-R {params.GRCh38} \
		-V {input.vcf} \
		-O {output.vcf} \
		--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
		--filter-name "my_filter"
		'''

rule cnvkit:
	input:
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		cns=RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.cns",
		bed=RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.bed"
	params:
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"]
	shell:
		'''
		# init bash
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		#source activate
		conda activate /datapool/yanzeqin/software/miniforge3/envs/cnvkit_old
		/datapool/yanzeqin/software/miniforge3/envs/cnvkit_old/bin/cnvkit.py batch {input.tumor}    --normal {input.normal} --fasta /datapool/yanzeqin/database/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa --annotate /datapool/yanzeqin/database/GRCh38_GENCODE/refFlat.txt -m wgs --output-dir {params.outputDir}
		cat {output.cns} |grep -vw chromosome |grep chr > {output.bed}
		'''

rule cnvkit_seg:
	input:
		cns=RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.cns"
	output:
		seg1=RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.seg",
		seg2=RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.filted.seg"
	shell:
		'''
		/datapool/yanzeqin/software/miniforge3/envs/cnvkit_old/bin/cnvkit.py export seg {input.cns} -o {output.seg1}
		cat {output.seg1}|grep chr > {output.seg2}
		'''

rule cnvkit_GISTIC:
	input:expand([RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.filted.seg"], sample=SAMPLES.keys())
	output:
		primaryseg= RESULT_DIR + "/" +config["output"]["cnvkit_GISTIC"] + "/primary_tumor.seg",
		result= RESULT_DIR + "/" +config["output"]["cnvkit_GISTIC"] + "/all_lesions.conf_90.txt"
	params:
		gisticdir= RESULT_DIR + "/" +config["output"]["cnvkit_GISTIC"]
	shell:
		'''
		mkdir -p {params.gisticdir}
		cat {input} |grep -vw ID > {output.primaryseg}
		/share/Data01/pengguoyu/bin/gistic2 -b {params.gisticdir} -seg {output.primaryseg} -refgene /share/Data01/pengguoyu/App/gistic2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat -broad 1 -brlen 0.5 -conf 0.95 -maxseg 10000 -qvt 0.1 -ta 0.21 -td 0.25 -genegistic 1
                '''


rule AmpliconArchitect1:
	input:
		bed=RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.bed",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		bed=RESULT_DIR + "/{sample}/" + config["output"]["AmpliconArchitect"] + "/{sample}.bed",
	params:
		out1=RESULT_DIR + "/{sample}/" + config["output"]["AmpliconArchitect"],
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["AmpliconArchitect"] + "/{sample}"
	shell:
		'''
		mkdir -p {params.out1}
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/yanzeqin/software/anaconda/envs/ampliconsuite
		python3 /datapool/yanzeqin/software/AmpliconArchitect/src/amplified_intervals.py --bed {input.bed} --bam {input.tumor} --out {params.outputDir} --ref GRCh38
		'''


rule AmpliconArchitect2:
	input:
		bed=RESULT_DIR + "/{sample}/" + config["output"]["AmpliconArchitect"] + "/{sample}.bed",
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		summry=RESULT_DIR + "/{sample}/" + config["output"]["AmpliconArchitect2"] + "/{sample}_summary.txt"
	params:
		out1=RESULT_DIR + "/{sample}/" + config["output"]["AmpliconArchitect2"],
		outputDir= RESULT_DIR + "/{sample}/" + config["output"]["AmpliconArchitect2"] + "/{sample}"
	shell:
		'''
		mkdir -p {params.out1}
		chmod a+r {params.out1}
		/datapool/yanzeqin/software/AA/AmpliconArchitect/docker/run_aa_docker.sh --bam {input.tumor} --bed {input.bed} --out {params.outputDir} --ref GRCh38 --downsample -1	
		#python3 /datapool/yanzeqin/software/AmpliconArchitect/src/AmpliconArchitect.py --bam {input.tumor} --bed {input.bed} --out {params.outputDir} --ref GRCh38 --downsample  -1 
		'''

rule AA_cnvkit:
	input:
		cns=RESULT_DIR + "/{sample}/" + config["output"]["cnvkit"] + "/{sample}_"+T+"_BQSR.call.cns",
		bam=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		out=RESULT_DIR + "/{sample}/" + config["output"]["AA_cnvkit"] + "/{sample}_AA_CNV_SEEDS.bed"
	params:
		id="{sample}",
		outdir= RESULT_DIR + "/{sample}/" + config["output"]["AA_cnvkit"]
	shell:
		'''
		mkdir -p {params.outdir}
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/yanzeqin/software/anaconda/envs/ampliconsuite
		python3 /datapool/yanzeqin/software/AmpliconSuite-pipeline/AmpliconSuite-pipeline.py -s {params.id} -t 4 --cnv_bed {input.cns} --bam {input.bam} --ref GRCh38 --output_directory {params.outdir}
		touch {output.out}
		'''



rule karken2_fastp:
	input:
		rawfq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
		rawfq2 = lambda wildcards: SAMPLES[wildcards.sample][1]
	output:
		read1 =  temp(RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"_karken_fastp1.fq.gz"),
		read2 =  temp(RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"_karken_fastp2.fq.gz"),
		fastp_html = RESULT_DIR + "/{sample}/" + config["output"]["fastp"] +  "/{sample}_"+T+"karken.report.html",
		fastp_json = RESULT_DIR + "/{sample}/" + config["output"]["fastp"] +  "/{sample}."+T+"karken.fastp.json"
	params:
		v1 = "{sample}_"+T,
	log: RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"karkenfastp.log"
	shell:
		r'''
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/fastp \
		-i {input.rawfq1} -o {output.read1} -I {input.rawfq2} -O {output.read2} \
		--json {output.fastp_json} --html {output.fastp_html} \
		--report_title {params.v1} \
		--detect_adapter_for_pe --compression 6 --cut_front --cut_tail   2> {log}
		'''

rule karken:
	input:
		read1 =  RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"_karken_fastp1.fq.gz",
		read2 =  RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_"+T+"_karken_fastp2.fq.gz"
	output:
		kra=RESULT_DIR + "/{sample}/" + config["output"]["kraken"] + "/{sample}_"+T+".kraken",
		kreport=RESULT_DIR + "/{sample}/" + config["output"]["kraken"] + "/{sample}_"+T+".kreport",
		bra=RESULT_DIR + "/{sample}/" + config["output"]["kraken"] + "/{sample}_"+T+".bracken"
	params:
		dir=RESULT_DIR + "/{sample}/" + config["output"]["kraken"]
	shell:
		'''
		mkdir -p {params.dir}
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/yanzeqin/software/anaconda/envs/kraken2
		kraken2 --paired --gzip-compressed --use-names --threads 8 --output {output.kra} --report {output.kreport} --db /datapool/zhuguanghui/multi.omics/stLFR/kraken2db   \
		{input.read1} {input.read2}
		bracken -d /datapool/zhuguanghui/multi.omics/stLFR/kraken2db -i {output.kreport} -o {output.bra} -t 8
		'''

###################################################################
##gatk的soamtic cnv流程
########
rule PreprocessIntervals:
	input:
		ref="/datapool/yanzeqin/database/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa"	
	output:
		inter="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Result/gatk/preprocessed_1000.interval_list"
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" PreprocessIntervals \
		-R {input.ref} \
		--bin-length 1000 \
		--padding 0 \
		--interval-merging-rule OVERLAPPING_ONLY \
		-O {output.inter}
		'''
rule AnnotateIntervals:
	input:
		ref="/datapool/yanzeqin/database/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa",
		inter="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Result/gatk/preprocessed_1000.interval_list"
	output:
		annointer="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Result/gatk/targets.preprocessed.1000.annotated.tsv"
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" AnnotateIntervals -L {input.inter} -R {input.ref} -imr OVERLAPPING_ONLY -O {output.annointer}
		'''

rule CollectFragmentCounts_T:
	input:
		inter="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Result/gatk/preprocessed_1000.interval_list",
		bam=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		hdf5=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_counts.hdf5"
	priority: 50
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" CollectReadCounts \
		-I {input.bam} \
		-L {input.inter} \
		--interval-merging-rule OVERLAPPING_ONLY \
		-O {output.hdf5}
		'''

rule CollectFragmentCounts_N:
	input:
		inter="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Result/gatk/preprocessed_1000.interval_list",
		bam=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam"
	output:
		hdf5=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+N+"_counts.hdf5"
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" CollectReadCounts \
		-I {input.bam} \
		-L {input.inter} \
		--interval-merging-rule OVERLAPPING_ONLY \
		-O {output.hdf5}
		'''

rule CreateReadCountPanelOfNormals:
	input:
		expand([RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+N+"_counts.hdf5"],sample=SAMPLES.keys()),
		annointer="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Result/gatk/targets.preprocessed.1000.annotated.tsv"
	output:
		pon="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Result/gatk/cnvponC.pon.hdf5"
	params:
		V=" -I ".join(expand([RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+N+"_counts.hdf5"],sample=SAMPLES.keys())),
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx6500m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" CreateReadCountPanelOfNormals \
		-I {params.V} \
		--annotated-intervals {input.annointer} \
		--minimum-interval-median-percentile 5.0 \
		--output {output.pon}
		'''


rule DenoiseReadCounts:
	input:
		hdf5=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_counts.hdf5",
		pon="/datapool/zhuguanghui/multi.omics/yanzeqin/FFPE/FFFPE2/Result/gatk/cnvponC.pon.hdf5"
	output:
		st=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_clean.standardizedCR.tsv",
		de=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_clean.denoisedCR.tsv"
	priority: 50
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" DenoiseReadCounts \
		-I {input.hdf5} \
		--count-panel-of-normals {input.pon} \
		--standardized-copy-ratios {output.st} \
		--denoised-copy-ratios {output.de}
		'''


rule ModelSegments:
	input:
		de=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_clean.denoisedCR.tsv"
	output:
		seg=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_clean.cr.seg"
	params:
		name="{sample}_"+T+"_clean",
		dir=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] 
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" ModelSegments \
		--denoised-copy-ratios {input.de} \
		--output {params.dir} \
		--output-prefix {params.name}
		ls {output.seg}
		'''

rule CallCopyRatioSegments:
	input:
		seg=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_clean.cr.seg"
	output:
		call=RESULT_DIR + "/{sample}/" + config["output"]["gatk_cnv"] + "/{sample}_"+T+"_clean.called.seg"
	shell:
		'''
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx64G -Xms16G" CallCopyRatioSegments \
		--input {input.seg} \
		--output {output.call}
		'''

rule frag:
	input:
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam"
	output:
		frag=RESULT_DIR + "/{sample}/" + config["output"]["JaBbA"] + "/cov.corrected.bw"
	params:
		outdir=RESULT_DIR + "/{sample}/" + config["output"]["JaBbA"]
	shell:
		'''
		mkdir -p {params.outdir}
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/zhuguanghui/multi.omics/yzqsoftware/miniconda/envs/JaBbA
		frag -b {input.tumor} -d /datapool/yanzeqin/database/frag/ -w 200 -o {params.outdir}
		'''

rule JaBbA:
	input:
		result= RESULT_DIR + "/{sample}/" + config["output"]["svaba"] + "/{sample}.svaba.somatic.sv.vcf",
		frag=RESULT_DIR + "/{sample}/" + config["output"]["JaBbA"] + "/cov.corrected.bw"
	output:
		jabba=RESULT_DIR + "/{sample}/" + config["output"]["JaBbA"] + "/jabba.simple.cnv.vcf"
	params:
		outdir=RESULT_DIR + "/{sample}/" + config["output"]["JaBbA"]
	priority: 50
	shell:
		'''	
		#conda init bash
		#source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		#conda activate /datapool/zhuguanghui/multi.omics/yzqsoftware/miniconda/envs/JaBbA
		#CPLEX_DIR=/opt/ibm/ILOG/CPLEX_Studio2211/
		jba {input.result} {input.frag} -o {params.outdir}
		#jba {input.result} {input.frag} -o {params.outdir}
		'''

rule smoove:
	input:
		tumor=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+T+"_BQSR.bam",
		normal=RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_"+N+"_BQSR.bam"
	output:
		vcf=RESULT_DIR + "/{sample}/" + config["output"]["smoove"] + "/{sample}-smoove.genotyped.vcf.gz"
	params:
		name =  "{sample}",
		outdir=RESULT_DIR + "/{sample}/" + config["output"]["smoove"]

	shell:
		'''
		conda init bash
		source /share/Data01/pengguoyu/App/anaconda3/etc/profile.d/conda.sh
		conda activate /datapool/zhuguanghui/multi.omics/yzqsoftware/miniconda/envs/smoove
		smoove call -x --name {params.name} --fasta /share/database/openData/GRCh38_GENCODE/GRCh38.primary_assembly.genome.fa -p 4 --genotype --outdir {params.outdir} {input.normal} {input.tumor}
		'''

