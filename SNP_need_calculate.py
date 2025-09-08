# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 16:48:23 2018

@author: Steven Han
"""


        
from itertools import islice
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-ld',help='eQTL_LD',metavar='')
parser.add_argument('-can',help='candidate_sQTL',metavar='')
parser.add_argument('-hmp',help='allforje.hmp.txt',metavar='')
parser.add_argument('-o',help='outputdir',metavar='')
args=parser.parse_args()
outdir=args.o
header=open('header','r').readline().strip().split()[11:]
transloc=set()
for i in open(args.ld,'r'):
    i=i.strip()
    if i.startswith('Zm00001eb'):
        i=i.split()
        i[1]=i[1].split('-')
        key1=i[0]+'-'+i[1][0]
        key2=i[0]+'-'+i[1][1]
        transloc.add(key1)
        transloc.add(key2)
qtl=[]
for i in open(args.can,'r'):
    i=i.strip().split()
    key=i[6]+'-'+i[0]+':'+i[9]
    if key in transloc:
        qtl.append([key,i[0],i[1],i[2]])
for i in range(len(qtl)):
    name=qtl[i][0].split('-')[0]+'C'+qtl[i][0].split('-')[1].split(':')[0]+'L'+qtl[i][0].split('-')[1].split(':')[1]
    out=open(outdir+name,'w')
    out.write('\t')
    out.write('\t'.join(header))
    out.write('\n')
    chr=qtl[i][1]
    start=qtl[i][2]
    end=qtl[i][3]
    for j in islice(open(args.hmp,'r'),1,None):
        j=j.strip().split()
        if j[2]==chr and int(j[3])>=int(start) and int(j[3])<=int(end):
            change={}
            allele=j[1].split('/')
            change[allele[0]]='0'
            change[allele[1]]='1'
            change['N']='NA'
            out.write(j[0])
            j=j[11:]
            for k in range(len(j)):
                out.write('\t'+change[j[k]])
            out.write('\n')
    out.close()
        
