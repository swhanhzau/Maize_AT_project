cat candidate5e9sQTLC candidate5e9sQTLD|cut -f6,7,8,10|sort|uniq >needJC
cut -f1 needJC|sort|uniq >jcsnp
cut -f2 needJC|sort|uniq >jctranscripts
rg -w -f jcsnp ~/zea_mays_gene_anno/sortV5SNP.hmp.txt >jcsnp_hmptemp
cat header jcsnp_hmptemp >jcsnp.hmp
rg -w -f jctranscripts matrix_wi_na_WW >jctranscripts_wwratiotemp
rg -w -f jctranscripts matrix_wi_na_DS >jctranscripts_dsratiotemp
cat ashead jctranscripts_wwratiotemp >jctranscripts_wwratio
cat ashead jctranscripts_dsratiotemp >jctranscripts_dsratio
python compute_r2.py --transcript jctranscripts_wwratio  --hmp jcsnp.hmp --snp_list needJC  --output JEresult_C
python compute_r2.py --transcript jctranscripts_dsratio  --hmp jcsnp.hmp --snp_list needJC  --output JEresult_D
