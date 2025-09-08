# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 15:56:55 2018

@author: Steven Han
"""
#usage: python eQTLsnpcluster.py -i significant_snp -o candidate_sQTL
#inputfile format: cluster\tmarker\tchr\tposition\tpvalue
#outputfile format: chr\tstart\tend\tcount\tleadP\tmarker\tcluster
#cluster is the phenotype's name
#the inputfile is from gwas result
from operator import itemgetter
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-i',help='allthreaholdpoint',metavar='')
parser.add_argument('-o',help='outputfile',metavar='')
args=parser.parse_args()
out=open(args.o,'w')
snpfile=args.i
snp=[]
for i in open(snpfile,'r'):
    i=i.strip().split()
    i[3]=int(i[3])
    i[4]=float(i[4])
    snp.append(i)
snp=sorted(snp, key=itemgetter(0,2,3))
start=end=-1
cluster=''
for i in range(len(snp)):
    chr_i=snp[i][2]
    loc_i=snp[i][3]-1
    if  (cluster==snp[i][0] and not( loc_i in range(start,end+1))) or (not cluster==snp[i][0]): #for tphenotypeis the cluster, we just care about snps in the same cluster
        leadp=snp[i][4]
        leadloc=loc_i+1
        marker=snp[i][1]
        cluster=snp[i][0]
        count=1
        start=loc_i
        end=loc_i+1
        for j in range(i+1,len(snp)):
            chr_j=snp[j][2]
            loc_j=snp[j][3]-1
            p=snp[j][4]
            if snp[i][0]==snp[j][0] and chr_i==chr_j and loc_j in range(start,end+10001):
                end=loc_j+1
                count+=1
                if p<leadp:
                    leadp=snp[j][4]
                    marker=snp[j][1]
                    leadloc=loc_j+1
                if j==len(snp)-1:
                    if count>=3:
                        print(chr_i,start,end,count,leadp,marker,cluster,chr_i,leadloc-1,leadloc,sep='\t',file=out)
            else:
                if count>=3:
                    print(chr_i,start,end,count,leadp,marker,cluster,chr_i,leadloc-1,leadloc,sep='\t',file=out)
                break

                
