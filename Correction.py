# Copyright © 2021, Shenghao Cao, Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

import sys
import os
import re
import argparse
from Bio import SeqIO
import numpy as np
import math


def CalcSeqLen(cigar):
    sections=re.split('([IMDSH])', cigar)
    length=0
    for k in range(len(sections)):
        if sections[k]=='M' or sections[k]=='D':
            length+=int(sections[k-1])
    return(length)

def recordScorePerBase(contigName, start, sections, seq, phred, trimL, trimR):
    moveRef=moveRead=0
    for k in range(len(sections)):
        if sections[k]=='M':
            L=0
            R=int(sections[k-1])
            if k==1 or k==3:
                L=trimL
            if k==len(sections)-2 or k==len(sections)-4:
                R=R-trimR
            for i in range(L,R):
                if seq[moveRead+i].upper() in scorePerBase[contigName][start+moveRef+i]:
                    scorePerBase[contigName][start+moveRef+i][seq[moveRead+i].upper()].append(1-10**(-(ord(phred[moveRead+i])-33)/10))
                if i > 0:
                    scorePerBase[contigName][start+moveRef+i]['insertion']['']+=1
            moveRef+=int(sections[k-1])
            moveRead+=int(sections[k-1])
        elif sections[k]=='D':
            for i in range(int(sections[k-1])):
                scorePerBase[contigName][start+moveRef+i][''].append(1)
                scorePerBase[contigName][start+moveRef+i]['insertion']['']+=1
            scorePerBase[contigName][start+moveRef+int(sections[k-1])]['insertion']['']+=1
            moveRef+=int(sections[k-1])
        elif sections[k]=='I':
            if seq[moveRead:moveRead+int(sections[k-1])] in scorePerBase[contigName][start+moveRef]['insertion']:
                scorePerBase[contigName][start+moveRef]['insertion'][seq[moveRead:moveRead+int(sections[k-1])]]+=1
            else:
                scorePerBase[contigName][start+moveRef]['insertion'][seq[moveRead:moveRead+int(sections[k-1])]]=1
            moveRead+=int(sections[k-1])
        elif sections[k]=='S':
            moveRead+=int(sections[k-1])

descript="This program corrects erroneous bases in the assembly.\n"
parser = argparse.ArgumentParser(description=descript)
parser.add_argument('-r1', required=True, help='Fastq file with forward paired reads')
parser.add_argument('-r2', required=True, help='Fastq file with reverse paired reads')
parser.add_argument('-t', type=int, default=1, help='Number of threads for parallelism')
parser.add_argument('-tr', type=int, default=2, help='Number of bases disregarded at alignment ends')
parser.add_argument('-i', default='complemented.fa', help='Input file')
parser.add_argument('-o', default='corrected.fa', help='Output file')
args = parser.parse_args()

print('\nMapping test reads to the draft assembly...\n')
os.system("bwa index %s 2>log"%args.i)
os.system("bwa mem -t %d %s %s > R1.sam 2>>log"%(args.t, args.i, args.r1))
os.system("bwa mem -t %d %s %s > R2.sam 2>>log"%(args.t, args.i, args.r2))
print('Correcting substitutions and indels...\n')

scorePerBase = {}
contigSeq={}

for seq_record in SeqIO.parse(args.i, 'fasta'):
    contigSeq[seq_record.id]=str(seq_record.seq)

with open('R1.sam') as fin_1:
    with open('R2.sam') as fin_2:
        record1_backup=['']*15;record2_backup=['']*15
        record1_num=0;record2_num=0

        while(True):

            while(True):
                record1=fin_1.readline().strip().split()

                if len(record1) == 3 and record1[0] == '@SQ':
                    contigName = record1[1].split(':')[1]
                    scorePerBase[contigName]={}

                    contigLen = int(record1[2].split(':')[1])
                    for i in range(contigLen):
                        scorePerBase[contigName][i+1]={'':[],'A':[],'G':[],'C':[],'T':[],'insertion':{'':0}}

                if len(record1) == 0:
                    break
                if not record1[0].startswith('@'):
                    if record1[0]!=record1_backup[0]:
                        break
                    else:
                        record1_num+=1

            while(True):
                record2=fin_2.readline().strip().split()

                if len(record2) == 0:
                    break
                if not record2[0].startswith('@'):
                    if record2[0]!=record2_backup[0]:
                        break
                    else:
                        record2_num+=1

            if record1_num==1:

                start = int(record1_backup[3])
                length=CalcSeqLen(record1_backup[5])
                if start < 300:
                    recordScorePerBase(record1_backup[2], start, re.split('([IMDSH])', record1_backup[5]), record1_backup[9], record1_backup[10], 0, args.tr)
                elif start+length > len(scorePerBase[record1_backup[2]])-300:
                    recordScorePerBase(record1_backup[2], start, re.split('([IMDSH])', record1_backup[5]), record1_backup[9], record1_backup[10], args.tr, 0)
                else:
                    recordScorePerBase(record1_backup[2], start, re.split('([IMDSH])', record1_backup[5]), record1_backup[9], record1_backup[10], args.tr, args.tr)

            if record2_num==1:

                start = int(record2_backup[3])
                length=CalcSeqLen(record2_backup[5])
                if start < 300:
                    recordScorePerBase(record2_backup[2], start, re.split('([IMDSH])', record2_backup[5]), record2_backup[9], record2_backup[10], 0, args.tr)
                elif start+length > len(scorePerBase[record2_backup[2]])-300:
                    recordScorePerBase(record2_backup[2], start, re.split('([IMDSH])', record2_backup[5]), record2_backup[9], record2_backup[10], args.tr, 0)
                else:
                    recordScorePerBase(record2_backup[2], start, re.split('([IMDSH])', record2_backup[5]), record2_backup[9], record2_backup[10], args.tr, args.tr)


            if len(record1)==0 and len(record2)==0:
                break
            else:
                record1_backup=record1;record2_backup=record2
                record1_num=1;record2_num=1

insertionList=[]
deletionList=[]
snpList=[]
trimList=[]

for contigName in scorePerBase:
    for base in scorePerBase[contigName]:

        insertion_dic=scorePerBase[contigName][base].pop('insertion')
        insertion=max(insertion_dic, key=insertion_dic.get).upper()

        A=G=C=T=D=0
        for prob in scorePerBase[contigName][base]['A']:
            A+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        for prob in scorePerBase[contigName][base]['G']:
            G+=math.log(0.99995*prob);A+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        for prob in scorePerBase[contigName][base]['C']:
            C+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);A+=math.log(0.99995*(1-prob)/3);T+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        for prob in scorePerBase[contigName][base]['T']:
            T+=math.log(0.99995*prob);G+=math.log(0.99995*(1-prob)/3);C+=math.log(0.99995*(1-prob)/3);A+=math.log(0.99995*(1-prob)/3);D+=math.log(0.00005/4)
        for prob in scorePerBase[contigName][base]['']:
            D+=math.log(1-0.00005);A+=math.log(0.00005);G+=math.log(0.00005);C+=math.log(0.00005);T+=math.log(0.00005)

        PA=1/(1+np.exp((T-A)*((T-A)<=709))+np.exp((G-A)*((G-A)<=709))+np.exp((C-A)*((C-A)<=709))+np.exp((D-A)*((D-A)<=709)))*((T-A)<=709)*((G-A)<=709)*((C-A)<=709)*((D-A)<=709)
        PG=1/(1+np.exp((A-G)*((A-G)<=709))+np.exp((T-G)*((T-G)<=709))+np.exp((C-G)*((C-G)<=709))+np.exp((D-G)*((D-G)<=709)))*((A-G)<=709)*((T-G)<=709)*((C-G)<=709)*((D-G)<=709)
        PC=1/(1+np.exp((A-C)*((A-C)<=709))+np.exp((G-C)*((G-C)<=709))+np.exp((T-C)*((T-C)<=709))+np.exp((D-C)*((D-C)<=709)))*((A-C)<=709)*((G-C)<=709)*((T-C)<=709)*((D-C)<=709)
        PT=1/(1+np.exp((A-T)*((A-T)<=709))+np.exp((G-T)*((G-T)<=709))+np.exp((C-T)*((C-T)<=709))+np.exp((D-T)*((D-T)<=709)))*((A-T)<=709)*((G-T)<=709)*((C-T)<=709)*((D-T)<=709)
        PD=1/(1+np.exp((A-D)*((A-D)<=709))+np.exp((G-D)*((G-D)<=709))+np.exp((C-D)*((C-D)<=709))+np.exp((T-D)*((T-D)<=709)))*((A-D)<=709)*((G-D)<=709)*((C-D)<=709)*((T-D)<=709)

        if PA==max(PA,PG,PC,PT,PD):
            site='A'
        elif PG==max(PA,PG,PC,PT,PD):
            site='G'
        elif PC==max(PA,PG,PC,PT,PD):
            site='C'
        elif PT==max(PA,PG,PC,PT,PD):
            site='T'
        elif PD==max(PA,PG,PC,PT,PD):
            site=''

        if A==G==C==T==D==0:
            trimList.append(base)
            site=contigSeq[contigName][base-1]
            insertion=''

        if site=='':
            insertionList.append(base)
        elif site != contigSeq[contigName][base-1]:
            snpList.append(1)
        if insertion!='':
            deletionList.append(len(insertion))

        scorePerBase[contigName][base]=insertion+site

    trimNum=0
    base=1
    while base in trimList:
        trimNum+=1
        scorePerBase[contigName][base]=''
        base+=1
    base=len(scorePerBase[contigName])
    while base in trimList:
        trimNum+=1
        scorePerBase[contigName][base]=''
        base-=1
# insertionNum=len(insertionList)
# for i in range(insertionNum-1):
#     if insertionList[i+1]-insertionList[i]==1:
#         insertionNum-=1

print("Corrected %d substitutions, inserted %d bases, deleted %d bases, trimmed %d bases.\n"%(sum(snpList),sum(deletionList),len(insertionList),trimNum))

with open(args.o,'w') as fout:
    contigNum=0
    for contigName in scorePerBase:
        contigNum+=1
        fout.write('>contig'+str(contigNum)+'\n')
        seq=('').join([scorePerBase[contigName][base] for base in scorePerBase[contigName]])
        line = int(len(seq)/100)
        for k in range(line):
            fout.write(seq[100*k:100*(k+1)]+'\n')
        if len(seq) > line*100:
            fout.write(seq[line*100:]+'\n')

os.system("rm %s.* log R1.sam R2.sam"%args.i)