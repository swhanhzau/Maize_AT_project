# Differential_Expression_and_Splicing_Analysis.py
# ==============================================
# Part 1: AT ttest (Transcript Ratio Differential Analysis)
# ==============================================
import numpy as np
from scipy import stats
import pandas as pd

# 1.1 Perform t-test on transcript ratio data and filter valid data
with open('matrix_wi_na_WW', 'r') as f1, open('matrix_wi_na_DS', 'r') as f2, \
        open('as_ttest_results.txt', 'w') as f_out, \
        open('valid_wwas.txt', 'w') as f_ww, open('valid_dsas.txt', 'w') as f_ds:
    
    # Write column names to result file
    f_out.write('GeneID\tTranscriptID\tMean_WW\tMean_DS\tDiff\tLog2FC\tp-value\n')
    
    # Read and write column names for valid data files (assuming first line is sample names)
    header_ww = f1.readline().strip()
    header_ds = f2.readline().strip()
    f_ww.write(f'{header_ww}\n')
    f_ds.write(f'{header_ds}\n')

    # Process WW and DS data line by line (ensuring both files have same number of lines)
    for line1, line2 in zip(f1, f2):
        line1 = line1.strip()
        line2 = line2.strip()
        if not line1 or not line2:  # Skip empty lines
            continue
        
        # Split transcript ID and expression ratio data
        parts1 = line1.split('\t')
        parts2 = line2.split('\t')
        transcript_id1, *data1_str = parts1
        transcript_id2, *data2_str = parts2
        
        # Extract gene ID from transcript ID (assuming format "GeneID_TranscriptID")
        gene_id1 = transcript_id1.split('_')[0] if '_' in transcript_id1 else transcript_id1
        gene_id2 = transcript_id2.split('_')[0] if '_' in transcript_id2 else transcript_id2
        
        # Data type conversion: string → float, "NA" → np.nan
        def convert_data(data_str):
            return [np.nan if x.strip().upper() == 'NA' else float(x.strip()) for x in data_str]
        data1 = convert_data(data1_str)
        data2 = convert_data(data2_str)
        
        # Calculate valid data ratio (number of samples excluding NaN and values > 0 / total samples)
        total_samples = len(data1)  # Assuming WW and DS have same number of samples
        valid1 = sum(1 for x in data1 if not np.isnan(x) and x > 0)
        valid2 = sum(1 for x in data2 if not np.isnan(x) and x > 0)
        valid_ratio1 = valid1 / total_samples if total_samples > 0 else 0
        valid_ratio2 = valid2 / total_samples if total_samples > 0 else 0
        
        # Only perform t-test on transcripts with valid ratio ≥50% in both WW and DS
        if valid_ratio1 < 0.5 or valid_ratio2 < 0.5:
            continue
        
        # Calculate statistical metrics
        mean_ww = np.nanmean(data1)
        mean_ds = np.nanmean(data2)
        diff = mean_ds - mean_ww
        # Calculate Log2FC (add 0.01 to avoid log(0))
        log2fc = np.log2((mean_ds + 0.01) / (mean_ww + 0.01)) if (mean_ww + 0.01 > 0 and mean_ds + 0.01 > 0) else np.nan
        # Paired t-test (ignore NaN)
        t_stat, p_val = stats.ttest_rel(data1, data2, nan_policy='omit')
        
        # Write t-test results
        f_out.write(f'{gene_id1}\t{transcript_id1}\t{mean_ww:.3f}\t{mean_ds:.3f}\t{diff:.3f}\t{log2fc:.3f}\t{p_val:.3g}\n')
        
        # Write valid data to corresponding files
        f_ww.write(f'{transcript_id1}\t' + '\t'.join([str(x) if not np.isnan(x) else 'NA' for x in data1]) + '\n')
        f_ds.write(f'{transcript_id2}\t' + '\t'.join([str(x) if not np.isnan(x) else 'NA' for x in data2]) + '\n')

print("Part 1 (AT ttest) completed. Output files: as_ttest_results.txt, valid_wwas.txt, valid_dsas.txt")


# ==============================================
# Part 2: exp ttest (Gene Expression Differential Analysis)
# ==============================================
# Data conversion function (redefined for Part 2)
def convert_data(data_str):
    return [np.nan if x.strip().upper() == 'NA' else float(x.strip()) for x in data_str]

with open('expww', 'r') as f1, open('expds', 'r') as f2, \
        open('exp_ttest_results.txt', 'w') as f_out, \
        open('valid_wwexp.txt', 'w') as f_ww, open('valid_dsexp.txt', 'w') as f_ds:
    
    # Write column names to result file
    f_out.write('GeneID\tMean_WW\tMean_DS\tDiff\tLog2FC\tp-value\n')
    
    # Read and write column names for valid data files (assuming first line is sample names)
    header_ww = f1.readline().strip()
    header_ds = f2.readline().strip()
    f_ww.write(f'{header_ww}\n')
    f_ds.write(f'{header_ds}\n')

    # Process WW and DS gene expression data line by line
    for line1, line2 in zip(f1, f2):
        line1 = line1.strip()
        line2 = line2.strip()
        if not line1 or not line2:
            continue
        
        # Split gene ID and expression data
        parts1 = line1.split('\t')
        parts2 = line2.split('\t')
        gene_id1, *data1_str = parts1
        gene_id2, *data2_str = parts2
        
        # Data type conversion (same as Part 1)
        data1 = convert_data(data1_str)
        data2 = convert_data(data2_str)
        
        # Calculate valid data ratio (number of samples excluding NaN and values > 0 / total samples)
        total_samples = len(data1)
        valid1 = sum(1 for x in data1 if not np.isnan(x) and x > 0)
        valid2 = sum(1 for x in data2 if not np.isnan(x) and x > 0)
        valid_ratio1 = valid1 / total_samples if total_samples > 0 else 0
        valid_ratio2 = valid2 / total_samples if total_samples > 0 else 0
        
        # Only perform t-test on genes with valid ratio ≥50% in either WW or DS
        if valid_ratio1 < 0.5 and valid_ratio2 < 0.5:
            continue
        
        # Calculate statistical metrics
        mean_ww = np.nanmean(data1)
        mean_ds = np.nanmean(data2)
        diff = mean_ds - mean_ww
        log2fc = np.log2((mean_ds + 0.01) / (mean_ww + 0.01)) if (mean_ww + 0.01 > 0 and mean_ds + 0.01 > 0) else np.nan
        t_stat, p_val = stats.ttest_rel(data1, data2, nan_policy='omit')
        
        # Write t-test results
        f_out.write(f'{gene_id1}\t{mean_ww:.3f}\t{mean_ds:.3f}\t{diff:.3f}\t{log2fc:.3f}\t{p_val:.3g}\n')
        
        # Write valid data to corresponding files
        f_ww.write(f'{gene_id1}\t' + '\t'.join([str(x) if not np.isnan(x) else 'NA' for x in data1]) + '\n')
        f_ds.write(f'{gene_id2}\t' + '\t'.join([str(x) if not np.isnan(x) else 'NA' for x in data2]) + '\n')

print("Part 2 (exp ttest) completed. Output files: exp_ttest_results.txt, valid_wwexp.txt, valid_dsexp.txt")


# ==============================================
# Part 3: Merge exp and ratio results (Mark DEG/DATG)
# ==============================================
# 3.1 Read t-test results from previous parts
df_exp = pd.read_table('exp_ttest_results.txt', sep='\t', na_values=['NA'])
df_as = pd.read_table('as_ttest_results.txt', sep='\t', na_values=['NA'])

# 3.2 Define filtering criteria for DEG and DAT (p-value correction: 0.05/total genes/total transcripts)
DEG_P_THRESHOLD = 0.05 / 27348  # 27348 is total number of genes
DEG_LOG2FC_THRESHOLD = 0.584963  # Corresponds to FC>1.5
DAT_P_THRESHOLD = 0.05 / 43798   # 43798 is total number of transcripts

# 3.3 Mark DEG (Differentially Expressed Genes)
df_exp['DEG'] = (df_exp['p-value'] < DEG_P_THRESHOLD) & (abs(df_exp['Log2FC']) > DEG_LOG2FC_THRESHOLD)

# 3.4 Mark DAT (Differentially Alternative Transcripts) and extract gene IDs for DATG
df_as['DAT'] = (df_as['p-value'] < DAT_P_THRESHOLD) & (df_as['p-value'].notna())
dat_genes = df_as[df_as['DAT']]['GeneID'].unique()  # Genes with at least one DAT transcript are DATG

# 3.5 Mark DATG in gene expression data
df_exp['DATG'] = df_exp['GeneID'].isin(dat_genes)

# 3.6 Merge gene expression and transcript ratio data (based on GeneID)
merged_df = pd.merge(
    df_exp, 
    df_as, 
    on='GeneID', 
    how='outer',  # Keep all genes/transcripts
    suffixes=('_gene', '_transcript')  # Distinguish same column names
)

# 3.7 Rename columns for better readability
merged_df.rename(columns={
    'Mean_WW_gene': 'GeneExp_WW',
    'Mean_DS_gene': 'GeneExp_DS',
    'Diff_gene': 'GeneExp_Diff',
    'Log2FC_gene': 'GeneExp_Log2FC',
    'p-value_gene': 'GeneExp_pval',
    'TranscriptID': 'TranscriptID',
    'Mean_WW_transcript': 'TransRatio_WW',
    'Mean_DS_transcript': 'TransRatio_DS',
    'Diff_transcript': 'TransRatio_Diff',
    'Log2FC_transcript': 'TransRatio_Log2FC',
    'p-value_transcript': 'TransRatio_pval'
}, inplace=True)

# 3.8 Calculate actual transcript expression (gene expression × transcript ratio / 100)
merged_df['TransExp_WW'] = merged_df['GeneExp_WW'] * merged_df['TransRatio_WW'] / 100
merged_df['TransExp_DS'] = merged_df['GeneExp_DS'] * merged_df['TransRatio_DS'] / 100
# Calculate Log2FC for transcript expression
merged_df['TransExp_Log2FC'] = np.log2(
    (merged_df['TransExp_DS'] + 0.01) / (merged_df['TransExp_WW'] + 0.01)
)

# 3.9 Format numeric values (keep 5 significant digits, avoid scientific notation redundancy)
def format_numeric(x):
    if isinstance(x, (int, float)) and not pd.isna(x):
        return "{:.5g}".format(x)
    return x
merged_df = merged_df.applymap(format_numeric)

# 3.10 Save merged results
merged_df.to_csv('merged_deg_datg.txt', sep='\t', index=False, na_rep='NA')

print("Part 3 (Merge Results) completed. Output file: merged_deg_datg.txt")
print("All analysis steps finished successfully!")