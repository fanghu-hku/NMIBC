#!/bin/bash
# run_example.sh
# 批量运行ShatterSeek分析示例
# 需要修改SAMPLES_FILE和PROJECT_DIR变量

SAMPLES_FILE="samples_all.csv"  # 格式: sample_name \t tumor_bam \t normal_bam
PROJECT_DIR="/path/to/project"  # 项目根目录

# 获取已完成样本
ls ${PROJECT_DIR}/*/ShatterSeek/*.chromSummary.txt 2>/dev/null | awk -F "/" '{print $8}' > sample_finish

# 运行未完成的样本
cat ${SAMPLES_FILE} | grep -vwf sample_finish | while read l1 l2 l3
do
    SAMPLE=${l1}
    OUTDIR=${PROJECT_DIR}/${SAMPLE}/ShatterSeek
    mkdir -p ${OUTDIR}
    cd ${OUTDIR}

    # 输入文件路径
    ANNO_MANTA=${PROJECT_DIR}/${SAMPLE}/neoSV/${SAMPLE}_manta.anno.txt
    ANNO_SVABA=${PROJECT_DIR}/${SAMPLE}/neoSV/${SAMPLE}_svaba.anno.txt
    ANNO_CNV=${PROJECT_DIR}/${SAMPLE}/JaBbA/jabba.seg

    SV_FILE=${OUTDIR}/sv.csv
    CN_FILE=${OUTDIR}/cn.csv

    # 准备输入文件
    echo -e "chrom\tstart\tend\tCN" > ${CN_FILE}
    cat ${ANNO_CNV} | grep -v "track.name" | grep -v KI | grep -v GL | grep -vw Y | grep -wv M | \
        awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$6}' >> ${CN_FILE}

    echo -e "chrom1\tstart1\tchrom2\tend2\tSVtype\tstrand1\tstrand2" > ${SV_FILE}
    cat ${ANNO_SVABA} ${ANNO_MANTA} | grep -v chrom1 | grep -vw "None" | grep -vw Y | grep -wv M | \
        awk -F "\t" '{print $1"\t"$2"\t"$8"\t"$9"\t"$16"\t"$6"\t"$13}' >> ${SV_FILE}

    # 运行ShatterSeek
    Rscript ShatterSeek.R -s ${SV_FILE} -c ${CN_FILE} -n ${SAMPLE} -o ${OUTDIR}

    echo "${SAMPLE} done"
done
