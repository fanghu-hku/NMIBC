#!/bin/bash
# sample.sh
# 根据样本名称判断Tumor/Normal类型
# 用法: sh sample.sh <input_file> <output_file>

input_file=$1
output_file=$2

# 读取输入文件并处理每一行
cat $input_file | while IFS=$'\t' read -r name path
do
    if [[ $name == *"N"* ]]; then
        echo -e "Normal\t${name}\t${path}"
    elif [[ $name == *"T"* ]]; then
        echo -e "Tumor\t${name}\t${path}"
    else
        echo -e "None\t${name}\t${path}"
    fi
done > $output_file

# 查看输出结果
head $output_file
