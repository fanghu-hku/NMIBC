#!/usr/bin/env Rscript
# ShatterSeek.R
# Chromothripsis detection analysis
# Usage: Rscript ShatterSeek.R -s sv.csv -c cn.csv -n sample_name -o output_dir

suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-s","--sv"), help = "sv file"),
    make_option(c("-c","--cnv"), help = "cnv bed file"),
    make_option(c("-n","--name"), help = "sample name"),
    make_option(c("-o","--outdir"), default = '.', help = "output directory")
)

opt_parser <- OptionParser(option_list = option_list,
                           description='\n\tThis script used to analyse chromothripsis.
                                        \nwrited by Zeqin Yan.')
args <- parse_args(opt_parser)

library(ShatterSeek)

SV <- read.table(args$sv, header=T, sep="\t")
CN <- read.table(args$cnv, header=T, sep="\t")
sample <- args$name
outdir1 <- args$outdir

outdir2 <- paste(outdir1, sample, sep="/")

print("finish 1")

SV_data <- SVs(chrom1=as.character(SV$chrom1),
               pos1=as.numeric(SV$start1),
               chrom2=as.character(SV$chrom2),
               pos2=as.numeric(SV$end2),
               SVtype=as.character(SV$SVtype),
               strand1=as.character(SV$strand1),
               strand2=as.character(SV$strand2))

print("finish 2")

CN_data <- CNVsegs(chrom=as.character(CN$chrom),
                   start=CN$start,
                   end=CN$end,
                   total_cn=CN$CN)

start_time <- Sys.time()

chromothripsis <- shatterseek(SV.sample=SV_data, seg.sample=CN_data, genome="hg38")

save(chromothripsis, file=paste(outdir2, "chromothripsis.Rdata", sep="."))

end_time <- Sys.time()

print(paste0("Running time (s): ", round(end_time - start_time, digits=2)))

write.table(chromothripsis@chromSummary, file=paste(outdir2, "chromSummary.txt", sep="."),
            sep="\t", col.names=T, row.names=F, quote=F)

chromothripsis@chromSummary[!is.na(chromothripsis@chromSummary$start), ]

library(gridExtra)
library(cowplot)
library(ggplot2)

# Plotting
for (ch in as.character(chromothripsis@chromSummary[!is.na(chromothripsis@chromSummary$start), ]$chrom)) {
    plots_chr <- plot_chromothripsis(ShatterSeek_output=chromothripsis, chr=ch, sample_name=sample, genome="hg38")
    plot_chr <- arrangeGrob(plots_chr[[1]], plots_chr[[2]], plots_chr[[3]], plots_chr[[4]],
                            nrow=4, ncol=1, heights=c(0.2,.4,.4,.4))
    ggsave(plot_chr, file=paste(outdir2, ch, "pdf", sep="."))
}
