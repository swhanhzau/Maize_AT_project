Maize Alternative Splicing and Expression QTL Project
Project Overview
This repository contains analysis pipelines for studying alternative splicing and expression quantitative trait loci (eQTL/atQTL) in maize. The project includes differential gene expression analysis, co-expression network analysis, RNA-seq processing for ZmMBF1, and e/atQTL mapping.

Analysis Workflows
1. Differential Expression and Alternative Splicing Analysis
Script: Differential_Expression_and_AT_Analysis.py

This pipeline identifies differentially expressed genes (DEGs) and genes with differential alternative transcript usage (DATG) from RNA-seq data.

2. Weighted Gene Co-expression Network Analysis (WGCNA)
Script: wgcna.r

Performs co-expression network analysis to identify modules of highly correlated genes and their associations with sample traits.

3. ZmMBF1 RNA-seq Data Processing
Processes RNA-seq data for the ZmMBF1 gene, including alignment and identification of differentially alternatively spliced transcripts.

3.1 Index Building
bash
sh build_index.sh
3.2 Alignment with HISAT2
bash
sh submit_align.sh
3.3 DATG Identification with DEXSeq
bash
Rscript dexseq.R
4. Expression and Alternative Splicing QTL (e/atQTL) Analysis
Identifies significant SNPs associated with expression and alternative splicing variations.

4.1 Extract Significant SNPs
bash
python filter_significant_p.py
4.2 Generate Candidate QTLs
bash
python sQTLsnpcluster.py -i significant_snp_C -o candidate_sQTL_C
python hotspotinput.py -b transcript.bed -c candidate_sQTL_C -o hotspotinput_C
perl HotSpotPlot_ver1.pl hotspotinput_C > sQTL_Cheatmap.svg
4.3 Filtering Using Linkage Disequilibrium (LD)
bash
cut -f6 candidate_sQTL_C | sort | uniq > allLeadC
grep -f allLeadC -F ~/zea_mays_gene_anno/sortV5SNP.hmp.txt -w > allleadC.hmp.txt
cat header allleadC.hmp.txt > useForLC.hmp.txt
run_pipeline.pl -Xmx80g -fork1 -h useForLC.hmp.txt -ld -ldType All -export LDoutC
cut -f1,2,7,8,14 LDoutC.txt | '$5>0.1{print}' > LDoutC_treat_filter0.1 
rm useForLC.hmp.txt LDoutC.txt allleadC.hmp.txt LDoutC_treat allLeadC
python LDtreatment.py -LD LDoutC_treat_filter0.1 -ca candidate5e9sQTLC -o sQTLC_LD
4.4 Joint Effect Analysis
bash
cut -f2 all1e2D5e9 | sort | uniq > allforjeC
grep -f allforjeC -F ~/zea_mays_gene_anno/sortV5SNP.hmp.txt -w > allforjeCinfo
cat header allforjeCinfo > allforjeC.hmp.txt
mkdir JointEffectunpre_C
python SNP_need_calculate.py -ld sQTLC_LD -can candidate5e9sQTLC -hmp allforjeC.hmp.txt -o JointEffectunpre_C/
mkdir preparedforJE_C
python jointprepareTreat.py -m matrix_wi_na_WW -un JointEffectunpre_C/ -o preparedforJE_C/
mkdir JE_Rout_C
for i in `ls -A preparedforJE_C`; do Rscript jointEffect.R preparedforJE_C/$i > JE_Rout_C/$i.out; done
python treatJEoutput.py -d JE_Rout_C/ -o JEresult_C
python finalfile_sQTL.py -j JEresult_C -l sQTLC_LD -o final_all_sQTLC
Usage
Each analysis can be run independently following the commands listed above. Please ensure all dependencies are installed and input files are properly formatted before execution.

Dependencies
Python 3

R

HISAT2

DEXSeq R package

TASSEL (for LD analysis)

License
[Specify License]

Contact
[Your Name/Team Name]
[Contact Information]
