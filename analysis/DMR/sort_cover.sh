#!/bin/bash
# sort_cover.sh
# 处理甲基化数据：去除chr前缀，调整坐标，计算甲基化比例，排序
# 用法: sh sort_cover.sh <input_file> <output_file> <temp_file> <bedtools_sorted> <chr_filtered>

input_file=$1
output_file=$2
temp_file=$3
temp_f2=$4
temp_f3=$5

# bedtools sort -i ${input_file} > ${temp_f2}
cat ${temp_f2} | grep chr > ${temp_f3}

# 使用 awk 处理数据
# 修改第一列的chr前缀
# 第二列数值减 1
# 第四列转换为甲基化比例
# 只保留前四列并输出到文件

awk '{
    gsub(/^chr/, "", $1);
    $2 = $2 - 1;
    $4 = $4 / 100;
    print $1"\t"$2"\t"$3"\t"$4;
}' OFS="\t" ${temp_f3} > $temp_file

# 对处理后的数据进行排序并输出到最终文件
sort -k1,1V -k2,2n $temp_file > $output_file

echo "finish 1"

# 清理临时文件
rm $temp_file
