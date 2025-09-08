#!/bin/bash
#SBATCH --job-name=build_index
#SBATCH --output=/DATA/storage/swhan/v5genome/build_index.out
#SBATCH --error=/DATA/storage/swhan/v5genome/build_index.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=your_partition

module load hisat2  # 加载 HISAT2 模块

# 创建 HISAT2 索引的目录
index_dir="/DATA/storage/swhan/v5genome/hisat2_index"
mkdir -p $index_dir

# 创建索引
hisat2-build -p 8 /DATA/storage/swhan/v5genome/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa $index_dir/Zea_mays_hisat2_index

