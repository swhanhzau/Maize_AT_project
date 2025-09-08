import pandas as pd
import argparse

# 创建一个解析器
parser = argparse.ArgumentParser(description='Process QTL data.')

# 添加期望的命令行参数
parser.add_argument('-c', '--candidateqtl_file', required=True, help='Path to the candidate QTL file')
parser.add_argument('-m', '--mainqtl_file', required=True, help='Path to the main QTL file')
parser.add_argument('-t', '--transcript_file', required=True, help='Path to the transcript file')
parser.add_argument('-o', '--output_file', default='merged_file.txt', help='Path to the output file')

# 解析命令行参数
args = parser.parse_args()

# 定义列名
candidateqtl_columns = ['qtlchr', 'qtlstart', 'qtlend', 'snpnum', 'leadP', 'leadsnp', 'transcriptID', 'leadchr', 'leadstart', 'leadend']
mainqtl_columns = candidateqtl_columns  # mainQTL文件和candidateQTL文件格式相同

# 读取 candidateqtl 和 mainqtl 数据，并指定列名
candidateqtl = pd.read_csv(args.candidateqtl_file, sep='\t', header=None, names=candidateqtl_columns)
mainqtl = pd.read_csv(args.mainqtl_file, sep='\t', header=None, names=mainqtl_columns)

# 将 candidateqtl 和 mainqtl 的每一行转换为字符串，然后进行比较
candidateqtl['QTL_Type'] = candidateqtl.astype(str).agg(' '.join, axis=1).isin(mainqtl.astype(str).agg(' '.join, axis=1)).replace({True: 'mainQTL', False: 'satelliteQTL'})

# 读取转录本文件，并指定列名（假设是标准的 BED 格式）
transcript_columns = ['chrom', 'chromStart', 'chromEnd', 'transcriptID', 'dot', 'strand']
transcript = pd.read_csv(args.transcript_file, sep='\t', header=None, names=transcript_columns)

# 将 candidateqtl 文件和转录本文件根据转录本 ID（即 candidateqtl 的第七列和 transcript 的第四列）进行合并
merged = pd.merge(candidateqtl, transcript, left_on='transcriptID', right_on='transcriptID', how='inner')

# 在最后一列前添加一列，其内容为 1
merged.insert(len(merged.columns) - 1, 'New_Column', 1)

# 保存合并后的数据
merged.to_csv(args.output_file, sep='\t', index=False, header=False)

print(f"Results saved to {args.output_file}")

