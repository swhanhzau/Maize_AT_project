# 加载必要的库
library("DEXSeq")

# 设置文件路径
count_files <- Sys.glob("/DATA/storage/swhan/mbf1rnaseq/first/*_sorted.bam")
sample_names <- sub("_sorted.bam", "", basename(count_files))
conditions <- c("cri23dt1", "cri23wtww1", "oe1_dt1", "oe1_ww1", "oe2_dt1", "oe2_ww1", "cri23dt2", "cri23wtww2", "oe1_dt2", "oe1_ww2", "oe2_dt2", "oe2_ww2", "cri23wtdt1", "cri23ww1", "cri23wtdt2", "cri23ww2")  # 根据实际样本命名适当调整

# 创建样本数据框
sampleData <- data.frame(sampleName = sample_names, condition = conditions, row.names = sample_names)

# 读取计数数据
flattened_gff <- "/DATA/storage/swhan/v5genome/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.51.dexseq.gff"
dxd <- DEXSeqDataSetFromHTSeq(countFiles,
                              sampleData = sampleData,
                              design = ~ sample + exon + condition:exon,
                              flattenedfile = flattened_gff)

# 归一化、估计大小因子和离散度
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)

# 进行差异可变剪接分析
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

# 获取结果
results <- DEXSeqResults(dxd)

# 保存结果
write.csv(as.data.frame(results), file = "/DATA/storage/swhan/mbf1rnaseq/first/dexseq_results.csv")

# 可视化特定基因
# dexseq_plot(dxd, featureID = "some_gene_id", expressionID = sampleData$sampleName, colorBy = "condition", legend = TRUE)

