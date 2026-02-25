#!/usr/bin/env Rscript
# run_sequenza.R
# Run sequenza for tumor purity and ploidy analysis
# Usage: Rscript run_sequenza.R <sample_name> <input_file> <output_dir>

library(sequenza)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
    cat("Usage: Rscript run_sequenza.R <sample_name> <input_file> <output_dir>\n")
    quit(status=1)
}

sample_name <- args[1]
input_file <- args[2]
outdir <- args[3]

# Extract data
test <- sequenza.extract(input_file, verbose=TRUE,
                         chromosome.list=c(paste0("chr", c(1:22, "X", "Y"))))

# Fit model
CP <- sequenza.fit(test)

# Output results
sequenza.results(sequenza.extract=test, cp.table=CP,
                 sample.id=sample_name, out.dir=outdir)
