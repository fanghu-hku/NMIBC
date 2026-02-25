# NMIBC Multi-omics Pipeline

Non-Muscle Invasive Bladder Cancer (NMIBC) multi-omics analysis pipeline.

## Directory Structure

```
NMIBC/
├── pipeline/          # Snakemake workflows
│   ├── DNA.smk        # WGS: variant calling, CNV, SV
│   ├── RNA.smk        # RNA-seq: fusion, quantification
│   ├── WGBS.smk       # Methylation analysis
│   └── smallRNA.smk   # miRNA/piRNA quantification
├── analysis/          # Downstream analysis scripts
│   ├── DMR/           # Differential methylation
│   ├── ShatterSeek/   # Chromothripsis detection
│   └── purity_ploidy/ # Tumor purity/ploidy
├── config/            # Configuration templates
│   ├── config.yaml.example
│   └── samples.tsv.example
└── README.md
```

## Quick Start

1. Copy and edit configuration:
```bash
cp config/config.yaml.example config/config.yaml
cp config/samples.tsv.example config/samples.tsv
# Edit paths in config.yaml and add samples to samples.tsv
```

2. Run pipeline:
```bash
# Dry run
snakemake -s pipeline/DNA.smk --configfile config/config.yaml -n

# Execute
snakemake -s pipeline/DNA.smk --configfile config/config.yaml -j 8
```

## Pipelines

| Pipeline | Description | Input |
|----------|-------------|-------|
| DNA.smk | WGS analysis: QC, alignment, variant calling (Mutect2, Strelka), CNV (FREEC, CNVkit, GATK CNV), SV (Manta, SvABA) | Tumor/Normal FASTQ |
| RNA.smk | RNA-seq: QC, STAR-Fusion, Salmon, TRUST4 | FASTQ |
| WGBS.smk | Methylation: Bismark alignment and extraction | FASTQ |
| smallRNA.smk | Small RNA: miRNA/piRNA quantification | FASTQ |

## Requirements

- Snakemake >= 5.0
- Conda environments for each tool
- Reference genomes and databases (see config.yaml.example)

## Author

Zeqin Yan
