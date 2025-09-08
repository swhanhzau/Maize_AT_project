#!/bin/bash
#SBATCH --job-name=alignment_$1
#SBATCH --output=/DATA/storage/swhan/mbf1rnaseq/first/%x_%j.out
#SBATCH --error=/DATA/storage/swhan/mbf1rnaseq/first/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G


# 样本名和文件路径
sample=$1
fq1=/DATA/storage/swhan/mbf1rnaseq/first/${sample}.clean.R1.fq.gz
fq2=/DATA/storage/swhan/mbf1rnaseq/first/${sample}.clean.R2.fq.gz

# 索引路径
index_dir="/DATA/storage/swhan/v5genome/hisat2_index/Zea_mays_hisat2_index"

# 运行 HISAT2 并通过管道将输出传递给 samtools 进行排序和索引
hisat2 -p 8 --dta -x $index_dir -1 $fq1 -2 $fq2 | samtools sort -@ 8 -o /DATA/storage/swhan/mbf1rnaseq/first/${sample}_sorted.bam
samtools index /DATA/storage/swhan/mbf1rnaseq/first/${sample}_sorted.bam
