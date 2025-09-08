# -*- coding: utf-8 -*-
"""
Created on Sun Sep 23 16:29:07 2018

@author: Steven Han
"""
from operator import itemgetter
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('-b','--transcriptbed',help='transcript.bed',metavar='')
parser.add_argument('-c','--candidate',help='candidate_sQTL_D',metavar='')
parser.add_argument('-o',help='outputfile',metavar='')
args=parser.parse_args()
out=open(args.o,'w')
transcript=args.transcriptbed
candidate=args.candidate
d={}
l=[]
for i in open(transcript,'r'):
    i=i.strip().split()
    d[i[3]]=i[0:3]
for i in open(candidate,'r'):
    i=i.strip().split()
    x=[i[6],d[i[6]],i[0:5]]
    x=(x[0]+'\t'+'\t'.join(x[1])+'\t'+'\t'.join(x[2])).split()
    x[-1]=float(x[-1])
    l.append(x)
l=sorted(l, key=itemgetter(8),reverse=True)
for i in l:
    i[-1]=str(i[-1])
    print('\t'.join(i),file=out)
