# Copyright © 2021, Shenghao Cao & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

import sys
import os
import re
import argparse
import numpy as np
import math
from Bio import SeqIO
from statsmodels.robust.scale import mad
#Chimeric reads are removed from primal reads when computing mapping rate

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
                # if i > 0:
                    # scorePerBase[contigName][start+moveRef+i]['insertion']['']+=1
            moveRef+=int(sections[k-1])
            moveRead+=int(sections[k-1])
        elif sections[k]=='D':
            for i in range(int(sections[k-1])):
                scorePerBase[contigName][start+moveRef+i][''].append(1)
                # scorePerBase[contigName][start+moveRef+i]['insertion']['']+=1
            # scorePerBase[contigName][start+moveRef+int(sections[k-1])]['insertion']['']+=1
            moveRef+=int(sections[k-1])
        elif sections[k]=='I':
            # scorePerBase[contigName][start+moveRef]['insertion']['I']+=1
            moveRead+=int(sections[k-1])
        elif sections[k]=='S':
            moveRead+=int(sections[k-1])


descript="This program evaluates the quality of assembly.\n"
parser = argparse.ArgumentParser(description=descript)
parser.add_argument('-r1', required=True, help='Fastq file with forward paired reads')
parser.add_argument('-r2', required=True, help='Fastq file with reverse paired reads')
parser.add_argument('-t', type=int, default=1, help='Number of threads for parallelism')
parser.add_argument('-i', default='assembly.fa', help='Input file')
# parser.add_argument('-o1', default='quality.perBase', help='Output file of quality value per base')
parser.add_argument('-o', default='summary.evaluation', help='Output file of evaluation statistics')
args = parser.parse_args()

print('\nMapping test reads to the polished assembly...\n')
os.system("bwa index %s 2>log"%args.i)
os.system("bwa mem -t %d %s %s > R1.sam 2>>log"%(args.t, args.i, args.r1))
os.system("bwa mem -t %d %s %s > R2.sam 2>>log"%(args.t, args.i, args.r2))

print('Computing the statistics of quality evaluation...\n')

contigSeq={}
for seq_record in SeqIO.parse(args.i, 'fasta'):
    contigSeq[seq_record.id]=str(seq_record.seq)

totalReads=mappedReads=multiReads=0

readCoverage = {}
scorePerBase = {}

with open('R1.sam') as fin_1:
    with open('R2.sam') as fin_2:
        record1_backup_list=[['']*15];record2_backup_list=[['']*15]
        record1_num=0;record2_num=0

        while(True):

            while(True):
                record1=fin_1.readline().strip().split()

                if len(record1) == 3 and record1[0] == '@SQ':
                    contigName = record1[1].split(':')[1]
                    scorePerBase[contigName]={}
                    readCoverage[contigName]={}

                    contigLen = int(record1[2].split(':')[1])
                    for i in range(contigLen):
                        readCoverage[contigName][i+1]=0
                        # scorePerBase[contigName][i+1]={'':[],'A':[],'G':[],'C':[],'T':[],'insertion':{'':0,'I':0}}
                        scorePerBase[contigName][i+1]={'':[],'A':[],'G':[],'C':[],'T':[]}

                if len(record1) == 0:
                    break
                if not record1[0].startswith('@'):
                    if record1[0]!=record1_backup_list[0][0]:
                        break
                    else:
                        if 'H' in record1[5]:
                            record1_num=-np.infty
                        else:
                            record1_num+=1
                            record1_backup_list.append(record1)

            while(True):
                record2=fin_2.readline().strip().split()

                if len(record2) == 0:
                    break
                if not record2[0].startswith('@'):
                    if record2[0]!=record2_backup_list[0][0]:
                        break
                    else:
                        if 'H' in record2[5]:
                            record2_num=-np.infty
                        else:
                            record2_num+=1
                            record2_backup_list.append(record2)

            if record1_num==1:
                record1_backup=record1_backup_list[0]
                if record1_backup[1]=='4':
                    totalReads+=1
                else:
                    totalReads+=1
                    mappedReads+=1
            elif record1_num>1:
                totalReads+=1
                mappedReads+=1
                multiReads+=1
            if record1_num>=1:
                for record1_backup in record1_backup_list:
                    start = int(record1_backup[3])
                    length=CalcSeqLen(record1_backup[5])
                    if start < 300:
                        recordScorePerBase(record1_backup[2], start, re.split('([IMDSH])', record1_backup[5]), record1_backup[9], record1_backup[10], 0, 2)
                        for i in range(start, start+length-2):
                            readCoverage[record1_backup[2]][i] += 1
                    elif start+length > len(scorePerBase[record1_backup[2]])-300:
                        recordScorePerBase(record1_backup[2], start, re.split('([IMDSH])', record1_backup[5]), record1_backup[9], record1_backup[10], 2, 0)
                        for i in range(start+2, start+length):
                            readCoverage[record1_backup[2]][i] += 1
                    else:
                        recordScorePerBase(record1_backup[2], start, re.split('([IMDSH])', record1_backup[5]), record1_backup[9], record1_backup[10], 2, 2)
                        for i in range(start+2, start+length-2):
                            readCoverage[record1_backup[2]][i] += 1

            if record2_num==1:
                record2_backup=record2_backup_list[0]
                if record2_backup[1]=='4':
                    totalReads+=1
                else:
                    totalReads+=1
                    mappedReads+=1
            elif record2_num>1:
                totalReads+=1
                mappedReads+=1
                multiReads+=1
            if record2_num>=1:
                for record2_backup in record2_backup_list:
                    start = int(record2_backup[3])
                    length=CalcSeqLen(record2_backup[5])
                    if start < 300:
                        recordScorePerBase(record2_backup[2], start, re.split('([IMDSH])', record2_backup[5]), record2_backup[9], record2_backup[10], 0, 2)
                        for i in range(start, start+length-2):
                            readCoverage[record2_backup[2]][i] += 1
                    elif start+length > len(scorePerBase[record2_backup[2]])-300:
                        recordScorePerBase(record2_backup[2], start, re.split('([IMDSH])', record2_backup[5]), record2_backup[9], record2_backup[10], 2, 0)
                        for i in range(start+2, start+length):
                            readCoverage[contigName][i] += 1
                    else:
                        recordScorePerBase(record2_backup[2], start, re.split('([IMDSH])', record2_backup[5]), record2_backup[9], record2_backup[10], 2, 2)
                        for i in range(start+2, start+length-2):
                            readCoverage[contigName][i] += 1

            if len(record1)==0 and len(record2)==0:
                break
            else:
                record1_backup_list=[record1];record2_backup_list=[record2]
                if 'H' in record1[5]:
                    record1_num=-np.infty
                else:
                    record1_num=1
                if 'H' in record2[5]:
                    record2_num=-np.infty
                else:
                    record2_num=1

mappingRate=round(mappedReads/totalReads,6)
multiRate=round(multiReads/totalReads,6)

readCoverageMedian = np.median([j for k in readCoverage.values() for j in k.values()])
readCoverageStd = mad([j for k in readCoverage.values() for j in k.values()])

#highReadCoverageNum = 0
lowReadCoverageNum = 0
totalBaseNum = 0
for contigName in readCoverage:
    #highReadCoverage=[]
    lowReadCoverage=[]
    for base in readCoverage[contigName]:
        #if readCoverage[contigName][base] > readCoverageMedian*2.8:
        #    highReadCoverage.append(base)
        if readCoverage[contigName][base] < readCoverageMedian-4*readCoverageStd:
            lowReadCoverage.append(base)

        # insertion_dic=scorePerBase[contigName][base].pop('insertion')
        # if sum(insertion_dic.values()) == 0 :
        #     aheadErrorRate = 1
        # else:
        #     I=insertion_dic['I']*math.log(1-0.00005)+insertion_dic['']*math.log(0.00005)
        #     D=insertion_dic['']*math.log(1-0.00005)+insertion_dic['I']*math.log(0.00005)
        #     aheadErrorRate = 1/(1+np.exp((D-I)*((D-I)<=708)))*((D-I)<=708)

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

        PA=1/(1+np.exp((T-A)*((T-A)<=708))+np.exp((G-A)*((G-A)<=708))+np.exp((C-A)*((C-A)<=708))+np.exp((D-A)*((D-A)<=708)))*((T-A)<=708)*((G-A)<=708)*((C-A)<=708)*((D-A)<=708)
        PG=1/(1+np.exp((A-G)*((A-G)<=708))+np.exp((T-G)*((T-G)<=708))+np.exp((C-G)*((C-G)<=708))+np.exp((D-G)*((D-G)<=708)))*((A-G)<=708)*((T-G)<=708)*((C-G)<=708)*((D-G)<=708)
        PC=1/(1+np.exp((A-C)*((A-C)<=708))+np.exp((G-C)*((G-C)<=708))+np.exp((T-C)*((T-C)<=708))+np.exp((D-C)*((D-C)<=708)))*((A-C)<=708)*((G-C)<=708)*((T-C)<=708)*((D-C)<=708)
        PT=1/(1+np.exp((A-T)*((A-T)<=708))+np.exp((G-T)*((G-T)<=708))+np.exp((C-T)*((C-T)<=708))+np.exp((D-T)*((D-T)<=708)))*((A-T)<=708)*((G-T)<=708)*((C-T)<=708)*((D-T)<=708)
        PD=1/(1+np.exp((A-D)*((A-D)<=708))+np.exp((G-D)*((G-D)<=708))+np.exp((C-D)*((C-D)<=708))+np.exp((T-D)*((T-D)<=708)))*((A-D)<=708)*((G-D)<=708)*((C-D)<=708)*((T-D)<=708)

        if contigSeq[contigName][base-1].upper()=='A':
            quality=PA
        elif contigSeq[contigName][base-1].upper()=='G':
            quality=PG
        elif contigSeq[contigName][base-1].upper()=='C':
            quality=PC
        elif contigSeq[contigName][base-1].upper()=='T':
            quality=PT
        else:
            quality=0

        # scorePerBase[contigName][base]={'quality': round(quality,6), 'aheadErrorRate': round(aheadErrorRate,6)}
        scorePerBase[contigName][base]=round(quality,6)

    # scorePerBase[contigName][1]['aheadErrorRate']=0

    totalBaseNum+=len(scorePerBase[contigName])

    k=1
    while k in lowReadCoverage:
        lowReadCoverage.remove(k)
        k+=1
    k=len(readCoverage[contigName])
    while k in lowReadCoverage:
        lowReadCoverage.remove(k)
        k-=1


    #highReadCoverageNum+=len(highReadCoverage)
    lowReadCoverageNum+=len(lowReadCoverage)

#highReadCoverageRate = round(highReadCoverageNum/totalBaseNum,6)
lowReadCoverageRate = round(lowReadCoverageNum/totalBaseNum,6)

qualityValue_list=[]

# aheadErrorRate_list=[]
# with open(args.o1,'w') as fout2:
#     table = "{0:<10}\t{1:<10}\t{2:<10}\t{3:<10}\n"
#     fout2.write(table.format("contigName", "base", "qualityValue", "aheadAccuracyRate"))
#     for contigName in scorePerBase:
#         for base in scorePerBase[contigName]:
#             qualityValue_list.append(scorePerBase[contigName][base]['quality'])
#             aheadErrorRate_list.append(scorePerBase[contigName][base]['aheadErrorRate'])
#             fout2.write(table.format(contigName, str(base), str(scorePerBase[contigName][base]['quality']), str(1-scorePerBase[contigName][base]['aheadErrorRate'])))

with open(args.o,'w') as fout:
    table = "{0:<35}\t{1:<10}\n"
    fout.write(table.format('Total contig number:', str(len(scorePerBase))))
    contigLenSeq=''
    for contigName in scorePerBase:
        contigLenSeq+=str(len(scorePerBase[contigName]))
        contigLenSeq+='\t'
        qualityValue_list+=scorePerBase[contigName].values()
    fout.write(table.format('Contig lengths:', contigLenSeq))
    fout.write(table.format('Mapping rate:', str(mappingRate*100)+'%'))
    fout.write(table.format('Multi-mapping rate:', str(multiRate*100)+'%'))
    # fout.write(table.format('High depth region:', str(highReadCoverageRate*100)+'%'))
    fout.write(table.format('Low depth fraction:', str(lowReadCoverageRate*100)+'%'))
    fout.write(table.format('Lowest quality value per base:', str(min(qualityValue_list))))
    # fout.write(table.format('Lowest accuracy ahead of each base:', str(1-max(aheadErrorRate_list))))

os.system("rm %s.* log R1.sam R2.sam"%args.i)
