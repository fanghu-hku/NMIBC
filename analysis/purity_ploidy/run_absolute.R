#!/usr/bin/env Rscript
# run_absolute.R
# Run ABSOLUTE for tumor purity and ploidy analysis
# Usage: Rscript run_absolute.R --seg sample.seg --name sample_name --outdir output_dir

suppressMessages(suppressWarnings(library(numDeriv)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(ABSOLUTE)))

option_list <- list(
    make_option("--seg", dest="seg.fn", help="segment file path"),
    make_option("--name", dest="sample.name", help="sample name"),
    make_option("--outdir", dest="results.dir", default=".", help="output directory"),
    make_option("--disease", dest="primary.disease", default="BLCA", help="primary disease type"),
    make_option("--min.ploidy", dest="min.ploidy", default="0.5", help="minimum ploidy"),
    make_option("--max.ploidy", dest="max.ploidy", default="8", help="maximum ploidy"),
    make_option("--max.sigma.h", dest="max.sigma.h", default="0.2", help="max sigma h"),
    make_option("--platform", dest="platform", default="SNP_6.0", help="platform type"),
    make_option("--copy.num.type", dest="copy.num.type", default="total", help="copy number type")
)

opt <- parse_args(OptionParser(option_list=option_list))

RunAbsolute(
    seg.dat.fn = opt$seg.fn,
    min.ploidy = as.numeric(opt$min.ploidy),
    max.ploidy = as.numeric(opt$max.ploidy),
    max.sigma.h = as.numeric(opt$max.sigma.h),
    platform = opt$platform,
    copy_num_type = opt$copy.num.type,
    sigma.p = 0,
    results.dir = opt$results.dir,
    primary.disease = opt$primary.disease,
    sample.name = opt$sample.name,
    max.as.seg.count = 1500,
    max.non.clonal = 1,
    max.neg.genome = 0.005
)
