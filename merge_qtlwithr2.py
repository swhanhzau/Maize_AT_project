import pandas as pd
import argparse

# 创建解析器
parser = argparse.ArgumentParser(description='Process QTL and JE data.')
parser.add_argument('je_file', type=str, help='The JE file to process.')
parser.add_argument('qtl_file', type=str, help='The QTL file to process.')
parser.add_argument('output_file', type=str, help='The file to output the merged data.')

# 解析参数
args = parser.parse_args()

# 读取 JE 文件
je_data = pd.read_csv(args.je_file, sep='\t', header=None, names=['transcriptID', 'SNPchr_loc', 'r2'])

# 修正 SNPchr_loc，确保 L 后面是整数
je_data['SNPchr_loc'] = je_data['SNPchr_loc'].str.replace(r'\.0', '', regex=True)  # 去掉 .0
je_data['SNPchr_loc'] = je_data['SNPchr_loc'].apply(lambda x: f"{x.split('L')[0]}L{int(float(x.split('L')[1]))}")  # 保证 L 后是整数

# 读取 QTL 数据
qtl_data = pd.read_csv(args.qtl_file, sep='\t', header=None, names=['qtlchr', 'qtlstart', 'qtlend', 'snpnum', 'leadP', 'leadsnp',
                                                              'transcriptID', 'SNPchr', 'leadstart', 'leadend', 'QTL_Type', 'transcriptchr',
                                                              'transcriptStart', 'transcriptEnd', 'name', 'score', 'strand'])

# 创建 SNPchr_loc 列，格式为 'CxLy'，与 JE 文件的格式一致
qtl_data['SNPchr_loc'] = 'C' + qtl_data['SNPchr'].astype(str) + 'L' + qtl_data['leadend'].astype(str)

# 合并 QTL 数据和 JE 数据
merged_data = pd.merge(qtl_data, je_data, on=['transcriptID', 'SNPchr_loc'], how='left')

# 判断 QTL 类型为 'cis' 或 'trans'
def classify_cis_trans(row):
    if str(row['qtlchr']) == str(row['transcriptchr']):
    # 如果 leadSNP 在转录本内，或者距离转录本 20kb 以内，就为 cis，否则为 trans
        if (row['leadend'] >= row['transcriptStart'] and row['leadend'] <= row['transcriptStart']) or \
           abs(row['leadstart'] - row['transcriptStart']) <= 20000 or \
           abs(row['leadend'] - row['transcriptEnd']) <= 20000:
            return 'cis'
        else:
            return 'trans'
    elif str(row['qtlchr']) != str(row['transcriptchr']):
        return 'diffchr'

# 应用函数来分类 'cis' 或 'trans'
merged_data['QTL_cis_trans'] = merged_data.apply(classify_cis_trans, axis=1)

# 输出合并后的数据到文件
merged_data.to_csv(args.output_file, sep='\t', index=False)

