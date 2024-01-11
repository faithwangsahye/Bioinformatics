#!/bin/bash
#组装+定量脚本

# 输入文件的目录
input_directory="/lustre/home/quq/data/hisat2"
# 输出文件的目录
output_directory="/lustre/home/quq/result2"
# GTF 文件路径
gtf_file="/lustre/home/quq/data/eRNA_hg38_2021.gff"

# 遍历目录中的所有 .bam 文件
for bam_file in "$input_directory"/*.bam; do
    # 获取不带路径和扩展名的文件名
    filename=$(basename -- "$bam_file")
    filename="${filename%.*}"

    # 构建输出文件名
    output_file="$output_directory/$filename.gtf"

    # 执行 stringtie 命令
    stringtie -p 3 -e -G "$gtf_file" -o "$output_file" -i "$bam_file"
done

python2 ~/stringtie-2.0.4.Linux_x86_64/prepDE.py \
-i /lustre/home/quq/data/sample_list.txt  \
-g gene_count_matrix.csv  \
-t transcript_count_matrix.csv
