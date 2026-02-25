#!/usr/bin/env Rscript
# plot_ShatterSeek.R
# 重新绘制ShatterSeek结果图
# 用法: Rscript plot_ShatterSeek.R -r rds_file -n sample_name -c chr -o output_dir --redir result_dir

suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-r","--rds"), help = "rds file"),
    make_option(c("-n","--name"), help = "sample name"),
    make_option(c("-c","--chr"), help = "chromosome"),
    make_option(c("-s","--start"), help = "start position"),
    make_option(c("-o","--outdir"), help = "output directory"),
    make_option("--redir", help = "result directory"),
    make_option(c("-e","--end"), help = "end position")
)

opt_parser <- OptionParser(option_list = option_list,
                           description='\n\tThis script used to plot chromothripsis results.
                                        \nwrited by Zeqin Yan.')
args <- parse_args(opt_parser)

sample <- args$name
outdir1 <- args$outdir
outdir2 <- paste(outdir1, sample, sep="/")
outdir3 <- args$redir
outdir4 <- paste(outdir3, sample, sep="/")
ch <- args$chr

print(ch)

library(ShatterSeek)

filep <- args$rds
load(filep)
print(filep)

library(gridExtra)
library(cowplot)
library(ggplot2)

plots_chr <- plot_chromothripsis(ShatterSeek_output=chromothripsis, chr=ch, sample_name=sample, genome="hg38")

plot_chr <- arrangeGrob(plots_chr[[1]], plots_chr[[2]], plots_chr[[3]], plots_chr[[4]],
                        nrow=4, ncol=1, heights=c(0.2,.4,.4,.4))

ggsave(plot_chr, file=paste(outdir4, ch, "_new.pdf", sep="."))
