#!/usr/bin/env Rscript
# run_sequenza.R
# 运行sequenza进行肿瘤纯度和倍性分析
# 用法: Rscript run_sequenza.R <sample_name> <input_file> <output_dir>

library(sequenza)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
    cat("Usage: Rscript run_sequenza.R <sample_name> <input_file> <output_dir>\n")
    quit(status=1)
}

sample_name <- args[1]
input_file <- args[2]
outdir <- args[3]

# 提取数据
test <- sequenza.extract(input_file, verbose=TRUE,
                         chromosome.list=c(paste0("chr", c(1:22, "X", "Y"))))

# 拟合模型
CP <- sequenza.fit(test)

# 输出结果
sequenza.results(sequenza.extract=test, cp.table=CP,
                 sample.id=sample_name, out.dir=outdir)
