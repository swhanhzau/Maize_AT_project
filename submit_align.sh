#!/bin/bash

# 定义样本数组
samples=("cri23dt1" "cri23wtww1" "oe1_dt1" "oe1_ww1" "oe2_dt1" "oe2_ww1" "cri23dt2" "cri23wtww2" "oe1_dt2" "oe1_ww2" "oe2_dt2" "oe2_ww2" "cri23wtdt1" "cri23ww1" "cri23wtdt2" "cri23ww2")

# 遍历所有样本
for sample in ${samples[@]}; do
    # 创建一个临时的 Slurm 脚本
    script_path="/DATA/storage/swhan/mbf1rnaseq/first/run_alignment_${sample}.sh"
    cat <<EOF >$script_path
#!/bin/bash
#SBATCH --job-name=alignment_$sample
#SBATCH --output=/DATA/storage/swhan/mbf1rnaseq/first/%x_%j.out
#SBATCH --error=/DATA/storage/swhan/mbf1rnaseq/first/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=Normal
# 样本文件路径
fq1=/DATA/storage/swhan/mbf1rnaseq/first/cleandata/${sample}.clean.R1.fq.gz
fq2=/DATA/storage/swhan/mbf1rnaseq/first/cleandata/${sample}.clean.R2.fq.gz

# 索引路径
index_dir="/DATA/storage/swhan/v5genome/hisat2_index/Zea_mays_hisat2_index"

# 运行 HISAT2 并通过管道将输出传递给 samtools 进行排序和索引
hisat2 -p 8 --dta -x $index_dir -1 \$fq1 -2 \$fq2 | samtools sort -@ 8 -o /DATA/storage/swhan/mbf1rnaseq/first/${sample}_sorted.bam
samtools index /DATA/storage/swhan/mbf1rnaseq/first/${sample}_sorted.bam
EOF

    # 提交任务
    sbatch $script_path
done

