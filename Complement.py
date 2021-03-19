# Copyright © 2021, Shenghao Cao & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import re


def rename_input(inp1, inp2):

    number=-1
    with open('testReads.fasta','w') as fout:
        for seq_record in SeqIO.parse(inp1, 'fastq'):
            number+=1
            fout.write('>'+str(number)+'\n')
            fout.write(str(seq_record.seq)+'\n')

        for seq_record in SeqIO.parse(inp2, 'fastq'):
            number+=1
            fout.write('>'+str(number)+'\n')
            fout.write(str(seq_record.seq)+'\n')


def joinContigs(consensus_list, minScore=20, trimLen=500):
    contigs = sorted(consensus_list, key=lambda k:len(k), reverse=True)
    contig_list = []
    index = list(range(len(contigs)))
    while len(index) > 0:
        ref = contigs[index.pop(0)]
        while True:
            maxScore = 0
            for j in index:
                query = contigs[j]
                alignments = aligner.align(ref[-trimLen:], query[:trimLen])
                alignment = alignments[0]
                if alignment.score>maxScore:
                    maxScore = alignment.score
                    maxQuery = query
                    maxAlignment = alignment
                    maxIndex = j
                alignments = aligner.align(ref[-trimLen:], str(Seq(query).reverse_complement())[:trimLen])
                alignment = alignments[0]
                if alignment.score>maxScore:
                    maxScore = alignment.score
                    maxQuery = str(Seq(query).reverse_complement())
                    maxAlignment = alignment
                    maxIndex = j
            if maxScore >= minScore:
                ref = ref[:-trimLen]+ref[-trimLen:][:maxAlignment.aligned[0][-1][-1]]+maxQuery[:trimLen][maxAlignment.aligned[1][-1][-1]:]+maxQuery[trimLen:]
                index.remove(maxIndex)
            else:
                break
        while True:
            maxScore = 0
            for j in index:
                query = contigs[j]
                alignments = aligner.align(query[-trimLen:],ref[:trimLen])
                alignment = alignments[0]
                if alignment.score>maxScore:
                    maxScore = alignment.score
                    maxQuery = query
                    maxAlignment = alignment
                    maxIndex = j
                alignments = aligner.align(str(Seq(query).reverse_complement())[-trimLen:], ref[:trimLen])
                alignment = alignments[0]
                if alignment.score>maxScore:
                    maxScore = alignment.score
                    maxQuery = str(Seq(query).reverse_complement())
                    maxAlignment = alignment
                    maxIndex = j
            if maxScore >= minScore:
                ref = maxQuery[:-trimLen]+maxQuery[-trimLen:][:maxAlignment.aligned[0][-1][-1]]+ref[:trimLen][maxAlignment.aligned[1][-1][-1]:]+ref[trimLen:]
                index.remove(maxIndex)
            else:
                break
        contig_list.append(ref)
    return(contig_list)

descript="This program joins separate contigs together.\n"
parser = argparse.ArgumentParser(description=descript)
parser.add_argument('-r1', required=True, help='fastq file with forward paired reads (required)')
parser.add_argument('-r2', required=True, help='fastq file with reverse paired reads (required)')
parser.add_argument('-i', default='draft.fa', help='input fasta file of draft assembly to be completed [default: draft.fa]')
parser.add_argument('-o', default='complemented.fa', help='fasta file to output complemented assembly to [default: complemented.fa]')
parser.add_argument('-t', type=int, default=1, help='number of threads for parallelism [default: 1]' )
args = parser.parse_args()

aligner = Align.PairwiseAligner(mode = 'local', match_score=1, mismatch_score=-2, open_gap_score=-5, extend_gap_score=-1)

rename_input(args.r1, args.r2)

os.system("bwa index %s 2>log"%args.i)
os.system("bwa mem -t 40 %s testReads.fasta > R.sam 2>>log"%args.i)

number = -1
reads=[]
contigLen={}

with open('complementaryReads.fa','w') as fout:
    for seq_record in SeqIO.parse(args.i, 'fasta'):
        number += 1
        fout.write('>'+str(number)+'\n')
        fout.write(str(seq_record.seq)+'\n')
        reads.append(str(seq_record.seq))
        contigLen[seq_record.id]=len(seq_record.seq)


    with open('R.sam') as fin:
        read_backup=''
        record_backup_list=[]
        
        for line in fin:
            record = line.strip().split()
            
            if not record[0].startswith('@'):
                if record[0]!=read_backup:
                    num=0
                    for record_backup in record_backup_list:
                        if record_backup[1] == '4':
                            num+=1
                        else:
                            sections=re.split('([IMDSH])', record_backup[5])
                            if (sections[1]=='S' and int(record_backup[3])<=5) or int(record_backup[3])>=contigLen[record_backup[2]]-len(record_backup[9]):
                                num+=1
                    if num==len(record_backup_list) and num in [1,2]:
                        number += 1
                        fout.write('>'+str(number)+'\n')
                        fout.write(record_backup[9]+'\n')
                        reads.append(record_backup[9])
                    record_backup_list=[record]
                    read_backup=record[0]
                else:
                    record_backup_list.append(record)

        num=0
        for record_backup in record_backup_list:
            if record_backup[1] == '4':
                num+=1
            else:
                sections=re.split('([IMDSH])', record_backup[5])
                if (sections[1]=='S' and int(record_backup[3])<=5) or int(record_backup[3])>=contigLen[record_backup[2]]-len(record_backup[9]):
                    num+=1
        if num==len(record_backup_list) and num in [1,2]:
            number += 1
            fout.write('>'+str(number)+'\n')
            fout.write(record_backup[9]+'\n')
            reads.append(record_backup[9])

os.system("rm %s.* log R.sam testReads.fasta"%args.i)
os.system("makeblastdb -in complementaryReads.fa -dbtype nucl -out db > /dev/null 2>&1")
os.system("rm -rf blastnResults")
os.system("mkdir blastnResults")
os.system("blastn -query complementaryReads.fa -db db -out blast.comp -evalue 0.01 -gapopen 5 -gapextend 2 -penalty -2 -reward 1 -word_size 10 -perc_identity 98 -outfmt 6")

graph = {}
selected = []
with open('blast.comp') as fin:
    for line in fin:
        record = line.strip().split()
        if int(record[0]) not in graph:
            graph[int(record[0])]={}
        if int(record[3])>=15 and record[0]!=record[1]:
            len1=len(reads[int(record[0])])
            len2=len(reads[int(record[1])])
            if (int(record[8])<int(record[9]) and min(int(record[6])-1, int(record[8])-1)<=5 and min(len1-int(record[7]), len2-int(record[9]))<=5) or (int(record[8])>int(record[9]) and min(int(record[6])-1, len2-int(record[8]))<=5 and min(len1-int(record[7]), int(record[9])-1)<=5):
                graph[int(record[0])][int(record[1])]={'alignmentLen':int(record[3]), 'start1':int(record[6]), 'end1':int(record[7]), 'start2':int(record[8]), 'end2':int(record[9])}

os.system("rm -rf db.* complementaryReads.fa blastnResults blast.comp")

finalList=[]
contigs=set(range(len(contigLen)))

while contigs!=set():
    root = list(contigs)[0]
    From=root
    selected.append(From)
    contigs.remove(From)
    finalSeq = reads[From]

    connect = 'tail'
    while True:
        candidate = {}
        if From in graph:
            for k in graph[From]:
                if k not in selected and ((connect == 'tail' and len(reads[From])-graph[From][k]['end1']<=5) or (connect == 'head' and graph[From][k]['start1']-1<=5)):
                    candidate[k]=graph[From][k]
        successor=list(contigs.intersection(set(candidate.keys())))
        if successor!=[]:
            To=successor[0]
            contigs.remove(To)
            selected.append(To)
        elif candidate!={}:
            To=min(candidate.keys(), key=lambda k:candidate[k]['alignmentLen'])
            selected.append(To)
        else:
            break
        if connect == 'tail' and graph[From][To]['start2']<graph[From][To]['end2']:
            finalSeq=finalSeq[:len(finalSeq)-(len(reads[From])-graph[From][To]['end1'])]+reads[To][graph[From][To]['end2']:]
            connect = 'tail'
        elif connect == 'head' and graph[From][To]['start2']>graph[From][To]['end2']:
            finalSeq=finalSeq[:len(finalSeq)-(graph[From][To]['start1']-1)]+reads[To][graph[From][To]['start2']:]
            connect = 'tail'
        elif connect == 'tail' and graph[From][To]['start2']>graph[From][To]['end2']:
            finalSeq=finalSeq[:len(finalSeq)-(len(reads[From])-graph[From][To]['end1'])]+str(Seq(reads[To]).reverse_complement())[len(reads[To])-graph[From][To]['end2']+1:]
            connect = 'head'
        elif connect == 'head' and graph[From][To]['start2']<graph[From][To]['end2']:
            finalSeq=finalSeq[:len(finalSeq)-(graph[From][To]['start1']-1)]+str(Seq(reads[To]).reverse_complement())[len(reads[To])-graph[From][To]['start2']+1:]
            connect = 'head'
        From=To

    From=root
    connect = 'head'
    while True:
        candidate = {}
        if From in graph:
            for k in graph[From]:
                if k not in selected and ((connect == 'tail' and len(reads[From])-graph[From][k]['end1']<=5) or (connect == 'head' and graph[From][k]['start1']-1<=5)):
                    candidate[k]=graph[From][k]
        successor=list(contigs.intersection(set(candidate.keys())))
        if successor!=[]:
            To=successor[0]
            contigs.remove(To)
            selected.append(To)
        elif candidate!={}:
            To=min(candidate.keys(), key=lambda k:candidate[k]['alignmentLen'])
            selected.append(To)
        else:
            break
        if connect == 'tail' and graph[From][To]['start2']<graph[From][To]['end2']:
            finalSeq=str(Seq(reads[To]).reverse_complement())[:len(reads[To])-graph[From][To]['end2']]+finalSeq[len(reads[From])-graph[From][To]['end1']:]
            connect = 'tail'
        elif connect == 'head' and graph[From][To]['start2']>graph[From][To]['end2']:
            finalSeq=str(Seq(reads[To]).reverse_complement())[:len(reads[To])-graph[From][To]['start2']]+finalSeq[graph[From][To]['start1']-1:]
            connect = 'tail'
        elif connect == 'tail' and graph[From][To]['start2']>graph[From][To]['end2']:
            finalSeq=reads[To][:graph[From][To]['end2']-1]+finalSeq[len(reads[From])-graph[From][To]['end1']:]
            connect = 'head'
        elif connect == 'head' and graph[From][To]['start2']<graph[From][To]['end2']:
            finalSeq=reads[To][:graph[From][To]['start2']-1]+finalSeq[graph[From][To]['start1']-1:]
            connect = 'head'
        From=To
    finalList.append(finalSeq)

contig_list = joinContigs(finalList)

with open(args.o,'w') as fout:
    contigNum=0
    for seq in contig_list:
        contigNum+=1
        fout.write('>contig'+str(contigNum)+'\n')
        line = int(len(seq)/100)
        for k in range(line):
            fout.write(seq[100*k:100*(k+1)]+'\n')
        if len(seq) > line*100:
            fout.write(seq[line*100:]+'\n')

table = "{0:<10}\t{1:<10}"
print(table.format('contigName', 'length'))
contigNum=0
for contig in contig_list:
    contigNum += 1
    print(table.format('contig%d'%contigNum, len(contig)))
print(table.format('\ntotal', sum([len(k) for k in contig_list])))
print('\n')