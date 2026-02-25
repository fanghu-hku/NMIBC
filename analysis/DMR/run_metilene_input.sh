#!/bin/bash
# run_metilene_input.sh
# 生成metilene输入文件
# 用法: sh run_metilene_input.sh <tumor_name_file> <normal_name_file> <input_dir>

TUMOR_FILE=$1
NORMAL_FILE=$2
INPUT_DIR=$3

# 初始化变量，不以逗号开头
Input1=""
Input2=""

# 读取癌症样本名称，并生成文件路径
while IFS= read -r l1 || [ -n "$l1" ]; do
    file1="${INPUT_DIR}/${l1}/FilterCoverage/Tumor_${l1}_Coverage_Filtered.txt"
    if [ -f "$file1" ]; then
        Input1="${Input1},${file1}"
    else
        echo "Warning: $file1 does not exist."
    fi
done < $TUMOR_FILE

# 读取正常样本名称，并生成文件路径
while IFS= read -r l1 || [ -n "$l1" ]; do
    file2="${INPUT_DIR}/${l1}/FilterCoverage/Normal_${l1}_Coverage_Filtered.txt"
    if [ -f "$file2" ]; then
        Input2="${Input2},${file2}"
    else
        echo "Warning: $file2 does not exist."
    fi
done < $NORMAL_FILE

# 移除字符串前面的逗号
Input1="${Input1#,}"
Input2="${Input2#,}"

# 运行metilene_input.pl生成输入文件
metilene_input.pl \
	--in1 ${Input1} \
	--in2 ${Input2} \
	--h1 Tumor --h2 Normal

# 运行后生成metilene_Tumor_Normal.input文件
