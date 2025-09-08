# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 11:21:11 2018

@author: Steven Han
"""

import argparse
import sys
parser=argparse.ArgumentParser()
parser.add_argument('-l',help='sQTLD_LD',metavar='')
parser.add_argument('-j',help='JEresult_D',metavar='')
parser.add_argument('-o',help='output',metavar='')
args=parser.parse_args()
sys.stdout=open(args.o,'w')
d={}
for i in open(args.j,'r'):
    i=i.strip().split()
    i[1]=i[1].split('L')[0][1:]+':'+str(int(float(i[1].split('L')[1])))
    key=i[0]+'-'+i[1]
    d[key]=float(i[2])
final=[]
need=[]
cal=[]
for i in open(args.l,'r'):
    i=i.strip().split()
    if len(i)==10:
        key=i[6]+'-'+i[7]+':'+i[9]
        if not d.get(key):
            final.append('\t'.join(i))
        else:
            need.append(i)
    elif len(i)==2:
        cal.append('-'.join(i))
for i in range(len(need)):
    if need[i] != 'deleted':
        clusteri=need[i][6]
        snpsumi=need[i][3]
        key1=need[i][7]+':'+need[i][9]
        for j in range(i+1,len(need)):
            if need[j] !='deleted' and need[j] !='deleted':
                clusterj=need[j][6]
                key2=need[j][7]+':'+need[j][9]
                snpsumj=need[j][3]
                if clusteri==clusterj:
                    if (clusteri+'-'+key1+'-'+key2 in cal) or (clusteri+'-'+key2+'-'+key1 in cal):
                        if d[clusteri+'-'+key1]>d[clusteri+'-'+key2]:
                            need[j]='deleted'
                        elif d[clusteri+'-'+key1]<d[clusteri+'-'+key2]:
                            need[i]='deleted'
                        elif d[clusteri+'-'+key1]==d[clusteri+'-'+key2]:
                            if snpsumi>=snpsumj:
                                need[j]='deleted'
                            else:
                                need[i]='deleted'
for i in range(len(need)):
    if need[i]!='deleted':
        final.append('\t'.join(need[i]))
for i in final:
    print(i)            
        
        
            
