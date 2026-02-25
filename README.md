# NMIBC Multi-omics Pipeline

Non-Muscle Invasive Bladder Cancer (NMIBC) multi-omics analysis pipeline.

## Directory Structure

```
NMIBC/
├── pipeline/                    # Snakemake workflows
│   ├── DNA.smk                  # WGS: variant calling, CNV, SV
│   ├── RNA.smk                  # RNA-seq: fusion, quantification
│   ├── WGBS.smk                 # Methylation analysis
│   └── smallRNA.smk             # miRNA/piRNA quantification
├── analysis/                    # Downstream analysis scripts
│   ├── DMR/                     # Differential methylation (metilene)
│   ├── ShatterSeek/             # Chromothripsis detection
│   └── purity_ploidy/           # Tumor purity/ploidy (ABSOLUTE, sequenza)
├── config/                      # Configuration file templates
│   ├── DNA.config.yaml.example
│   ├── RNA.config.yaml.example
│   ├── WGBS.config.yaml.example
│   ├── smallRNA.config.yaml.example
│   └── cluster.yaml.example
└── README.md
```

## Quick Start

1. Copy and edit configuration:
```bash
# For DNA pipeline
cp config/DNA.config.yaml.example config/DNA.config.yaml
cp config/DNA.samples.tsv.example config/DNA.samples.tsv

# Edit paths in config file and add samples
```

2. Run pipeline:
```bash
# Dry run
snakemake -s pipeline/DNA.smk --configfile config/DNA.config.yaml -n

# Execute
snakemake -s pipeline/DNA.smk --configfile config/DNA.config.yaml -j 8

# With cluster support
snakemake -s pipeline/DNA.smk --configfile config/DNA.config.yaml \
    --cluster "qsub -cwd -q all.q -l vf={cluster.mem}" \
    --cluster-config config/cluster.yaml -j 100
```

## Pipelines

| Pipeline | Description | Input | Config |
|----------|-------------|-------|--------|
| DNA.smk | WGS: QC, alignment, variant calling (Mutect2, Strelka), CNV (FREEC, CNVkit, GATK), SV (Manta, SvABA) | Tumor/Normal FASTQ | DNA.config.yaml |
| RNA.smk | RNA-seq: QC, STAR-Fusion, Salmon, TRUST4, Kraken2 | FASTQ | RNA.config.yaml |
| WGBS.smk | Methylation: Bismark alignment and extraction | FASTQ | WGBS.config.yaml |
| smallRNA.smk | Small RNA: miRNA/piRNA quantification | FASTQ (single-end) | smallRNA.config.yaml |

## Sample List Format

| Pipeline | Format |
|----------|--------|
| DNA | `sample_name\ttumor_fq1\ttumor_fq2\tnormal_fq1\tnormal_fq2` |
| RNA | `sample_name\tfq1\tfq2` |
| WGBS | `sample_name\tfq1\tfq2` |
| smallRNA | `sample_name\tfq` (single-end) |

## Requirements

- Snakemake >= 5.0
- Conda environments for each tool
- Reference genomes and databases (see config files)

## Author

Zeqin Yan
