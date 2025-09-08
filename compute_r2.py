# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 17:45:22 2024

@author: swhan
"""

import pandas as pd
from scipy.stats import linregress
import argparse

def load_and_prepare_transcript(transcript_file):
    """加载转录本比率文件并转置"""
    transcript_data = pd.read_csv(transcript_file, sep='\t')
    # 转置数据以便转录本为列，样本或条件为行
    transcript_data_transposed = transcript_data.set_index('transcript_id').transpose()
    return transcript_data_transposed

def load_and_prepare_hmp(hmp_file):
    """加载HMP文件，跳过第1到第10列，保留第0列和第11列以后的列"""
    with open(hmp_file, 'r') as file:
        columns = file.readline().strip().split('\t')
    cols_to_use = [0] + list(range(11, len(columns)))  # 第0列 + 第11列及以后的列
    hmp_data = pd.read_csv(hmp_file, sep='\t', usecols=cols_to_use, index_col=0)
    hmp_data_transposed = hmp_data.transpose()  # 转置数据
    return hmp_data_transposed

def convert_snp_to_binary(snp_series):
    """将SNP数据中两种碱基分别转换为0和1，删除包含N的行"""
    # 获取唯一值（排除N）
    unique_values = snp_series.unique()
    unique_values = [base for base in unique_values if base != 'N']
    
    # 确保只剩两种碱基
    if len(unique_values) == 2:
        base_mapping = {unique_values[0]: 0, unique_values[1]: 1}
        return snp_series.map(base_mapping)
    else:
        return None  # 如果不满足条件，返回None

def join_and_compute_r2(transcript_data, hmp_data, snp_list_file):
    """联立数据并计算R2"""
    snp_list = pd.read_csv(snp_list_file, sep='\t', names=['snp','transcript_id', 'chr', 'pos'])
    results = []

    for idx, row in snp_list.iterrows():
        transcript_id = row['transcript_id']
        snp = row['snp']
        chr_pos = f"C{row['chr']}L{row['pos']}"  # 构造 CxLy 格式

        if snp in hmp_data.columns and transcript_id in transcript_data.columns:
            # 提取对应的SNP数据和转录本数据
            snp_data = hmp_data[snp]
            transcript_expr = transcript_data[transcript_id]
            
            # 转换SNP数据为0和1，删除为N的行
            snp_data_binary = convert_snp_to_binary(snp_data)
            
            if snp_data_binary is not None:
                combined_data = pd.DataFrame({
                    'SNP': snp_data_binary,
                    'Transcript_Ratio': transcript_expr
                }).dropna()

                if not combined_data.empty:
                    # 计算R2
                    slope, intercept, r_value, p_value, std_err = linregress(
                        combined_data['SNP'], combined_data['Transcript_Ratio']
                    )
                    r2 = r_value ** 2  # R-squared

                    # 保存结果
                    results.append({
                        'Transcript_ID': transcript_id,
                        'SNP_Pos': chr_pos,
                        'R2': round(r2, 5)  # 保留5位小数
                    })

    # 返回结果
    return pd.DataFrame(results)

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="Compute R2 values for SNP-Transcript pairs")
    parser.add_argument('--transcript', required=True, help="Path to transcript ratio file")
    parser.add_argument('--hmp', required=True, help="Path to HMP file")
    parser.add_argument('--snp_list', required=True, help="Path to SNP list file")
    parser.add_argument('--output', required=True, help="Path to output file")
    args = parser.parse_args()

    # 加载并准备数据
    transcript_data = load_and_prepare_transcript(args.transcript)
    hmp_data = load_and_prepare_hmp(args.hmp)

    # 联立数据并计算R2
    final_data = join_and_compute_r2(transcript_data, hmp_data, args.snp_list)

    # 输出结果到文件
    final_data.to_csv(args.output, sep='\t', index=False, header=False)

if __name__ == '__main__':
    main()

