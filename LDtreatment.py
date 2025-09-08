# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 21:42:11 2018

@author: Steven Han
"""
from itertools import islice
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-LD',help='LDoutD_treat',metavar='')
parser.add_argument('-ca',help='candidate_sQTL_D',metavar='')
parser.add_argument('-o',help='eQTLD_LD',metavar='')
args=parser.parse_args()
outfile=open(args.o,'w')
ld=set()
qtl=[]
cluster=[]
key=[]
p=[]
need=set()
for i in islice(open(args.LD,'r'),1,None):
    i=i.strip().split()
    key1=i[0]+':'+i[1]
    key2=i[2]+':'+i[3]
    ld.add(key1+'-'+key2)
    ld.add(key2+'-'+key1)
for i in open(args.ca,'r'):
    i=i.strip()
    qtl.append(i)
    i=i.split()
    cluster.append(i[6])
    key.append(i[7]+':'+i[9])
    p.append(float(i[4]))
for i in range(len(qtl)-1):
    for j in range(i+1,len(qtl)):
        if cluster[i]==cluster[j]:
            if key[i]+'-'+key[j] in ld:
                if p[i]<p[j]:
                    qtl[j]='deleted'
                elif p[i]>p[j]:
                    qtl[i]='deleted'
                else:
                    need.add(cluster[i]+'__'+key[i]+'-'+key[j])
for i in need:
    print(i.replace('__','\t'),file=outfile)
for i in qtl:
    if i != 'deleted':
        print(i,file=outfile)
                    