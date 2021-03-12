# Copyright © 2021, Shenghao Cao & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

import os
import argparse
from Bio import SeqIO
import collections
import sys
import numpy as np
import random

def sampleReads(prefix, number):
    index = random.sample(range(int(len(total_reads_1)/4)), number)
    with open(prefix+"-1.fq",'w') as fout:
        for k in index:
            fout.write(total_reads_1[4*k])
            fout.write(total_reads_1[4*k+1])
            fout.write(total_reads_1[4*k+2])
            fout.write(total_reads_1[4*k+3])
    with open(prefix+"-2.fq",'w') as fout:
        for k in index:
            fout.write(total_reads_2[4*k])
            fout.write(total_reads_2[4*k+1])
            fout.write(total_reads_2[4*k+2])
            fout.write(total_reads_2[4*k+3])


descript="This program assembles genomes via robust regression and re-sampling techniques.\n"
parser = argparse.ArgumentParser(description=descript)
parser.add_argument('-t', type=int, default=1, help='number of threads for parallelism [default: 1]')
parser.add_argument('-s', type=int, default=10, help='total sampling times [default: 1]')
parser.add_argument('-asm', default="RegAssembler", help='alternative assemblers: RegAssembler, SPAdes [default: RegAssembler]')
parser.add_argument('-r1', default='filteredReads1.fq', help='fastq file with forward paired reads [default: filteredReads1.fq]')
parser.add_argument('-r2', default='filteredReads2.fq', help='fastq file with reverse paired reads [default: filteredReads2.fq]')
parser.add_argument('-n1', type=int, default=10000, help='number of training reads in each sampling for generating draft assembly [default: 10000]')
parser.add_argument('-n2', type=int, default=20000, help='number of test reads in each sampling for polishing and evaluating draft assembly [default: 20000]')
parser.add_argument('-thr', type=int, default=3, help='convergence threshold for modified IRLS algorithm to stop iteration [default: 3]')
parser.add_argument('-ho', type=int, default=2, help='admissible hanging-out length for each pair of overlapping reads [default: 2]')
parser.add_argument('-al', type=int, default=20, help='minimum alignment length for a successful overlap [default: 20]')
parser.add_argument('-nchi', action="store_true", help='this flag cancels chimera detection and removal')

args = parser.parse_args()

os.system("rm -rf log* out*")
os.system("rm -rf reads")
os.mkdir('reads')
if args.asm.lower() == "regassembler" or args.asm.lower() == "spades":
    os.system("rm -rf candidates")
    os.mkdir('candidates')
    d="candidates"
else:
    print('Please choose one of the following assemblers: RegAssembler, SPAdes.')
    sys.exit()

n1 = int(args.n1/2)
n2 = int(args.n2/2)
print('\nImporting read files ...\n')
readFile_1 = open(args.r1,'r')
readFile_2 = open(args.r2,'r')
total_reads_1 = readFile_1.readlines()
total_reads_2 = readFile_2.readlines()
readFile_1.close()
readFile_2.close()

for i in range(args.s):
    i+=1
    print('\n----------------------------------------------------------------------------')
    print('\nAssemble the %d/%d candidate genome.\n'%(i,args.s))
    n=0
    success = False

    while not success:
        if n<20:
            n+=1
        else:
            print('Failed to assemble a qualified genome,\nyou can change the sampling coverage or tune the parameters.')
            sys.exit()

        print('--------------------------------------')
        print('Produce the draft assembly.\n')
        print('Sampling training reads ...\n')
        sampleReads('reads/train%d'%i, n1)
        if args.asm.lower() == "regassembler":
            print('Running RegAssembler ...\n')
            if args.nchi:
                os.system("cd %s; python ../RegAssembler.py -r1 ../reads/train%d-1.fq -r2 ../reads/train%d-2.fq -nchi -t %d -thr %d -ho %d -al %d -o draft-%d.fa "%(d,i,i,args.t,args.thr,args.ho,args.al,i))
            else:
                os.system("cd %s; python ../RegAssembler.py -r1 ../reads/train%d-1.fq -r2 ../reads/train%d-2.fq -t %d -thr %d -ho %d -al %d -o draft-%d.fa "%(d,i,i,args.t,args.thr,args.ho,args.al,i))
        elif args.asm.lower() == "spades":
            print('Running SPAdes ...\n')
            os.system("cd %s; spades.py -1 ../reads/train%d-1.fq -2 ../reads/train%d-2.fq -o out-%d --isolate -t %d >> spades.log 2>&1"%(d,i,i,i,args.t))
            contigNum=0
            with open('%s/draft-%d.fa'%(d,i),'w') as fout:
                for seq_record in SeqIO.parse('%s/out-%d/contigs.fasta'%(d,i), 'fasta'):
                    if len(str(seq_record.seq))>500:
                        contigNum+=1
                        fout.write(">contig%d\n"%contigNum)
                        fout.write(str(seq_record.seq)+'\n')

        totalContigLen=0
        for seq_record in SeqIO.parse('%s/draft-%d.fa'%(d,i), 'fasta'):
            totalContigLen+=len(str(seq_record.seq))

        if totalContigLen>28400:
            m=0
            while m<20:
                m+=1

                print('--------------------------------------')
                print('Complete the draft assembly.\n')
                print('Sampling test reads ...\n')
                sampleReads('reads/test%d'%i, n2)

                print('Joining separate contigs ...\n')
                os.system("cd %s; python ../Complement.py -r1 ../reads/test%d-1.fq -r2 ../reads/test%d-2.fq -t %d -i draft-%d.fa -o complemented-%d.fa"%(d,i,i,args.t,i,i))
                maxContigLen=0
                maxContigSeq=''
                for seq_record in SeqIO.parse('%s/complemented-%d.fa'%(d,i), 'fasta'):
                    if len(str(seq_record.seq))>maxContigLen:
                        maxContigLen=len(str(seq_record.seq))
                        maxContigSeq=str(seq_record.seq)
                if maxContigLen>29500:
                    with open('%s/complete-%d.fa'%(d,i),'w') as fout:
                        fout.write(">contig1\n")
                        fout.write(maxContigSeq+'\n')

                    print('--------------------------------------')
                    print('Run the 1st round of error correction.\n')
                    os.system("cd %s; python ../Correction.py -r1 ../reads/test%d-1.fq -r2 ../reads/test%d-2.fq -i complete-%d.fa -o corrected-%d.fa -t %d -tr 15"%(d,i,i,i,i,args.t))

                    print('--------------------------------------')
                    print('Run the 2nd round of error correction.\n')
                    os.system("cd %s; python ../Correction.py -r1 ../reads/test%d-1.fq -r2 ../reads/test%d-2.fq -i corrected-%d.fa -o assembly-%d.fa -t %d -tr 2"%(d,i,i,i,i,args.t))

                    print('--------------------------------------')
                    print('Perform quality evaluation.\n')
                    os.system("cd %s; python ../Evaluation.py -r1 ../reads/test%d-1.fq -r2 ../reads/test%d-2.fq -i assembly-%d.fa -o summary.evaluation.%d"%(d,i,i,i,i))
                    os.system("cat %s/summary.evaluation.%d"%(d,i))

                    with open("%s/summary.evaluation.%d"%(d,i)) as fin:
                        for line in fin:
                            record=line.strip().split()
                            if record[0]=='Mapping':
                                mappingRate=float(record[-1].split("%")[0])*0.01
                            if record[0]=='Multi-mapping':
                                multiMappingRate=float(record[-1].split("%")[0])*0.01
                            if record[0]=='Low':
                                lowRate=float(record[-1].split("%")[0])*0.01
                            if record[1]=='quality':
                                lowestQ=float(record[-1])
                    if mappingRate>0.95 and multiMappingRate==0 and lowRate==0 and lowestQ>0.9:
                        success = True
                    break

print('\n\n----------------------------------------------------------------------------')
print('Determine the final assembly by mutiple sequence alignment.\n')
# Generate consensus
os.system("rm -rf msa")
os.mkdir('msa')


os.system('cd msa; makeblastdb -in ../%s/assembly-1.fa -dbtype nucl -out db > log'%d)
with open('msa/msaInput','w') as fout:
    for i in range(args.s):
        i+=1
        fout.write('>assembly%d\n'%i)
        os.system('cd msa; blastn -query ../%s/assembly-%d.fa -out blast.out-%d -db db -outfmt 6 -evalue 1e-5 -num_threads %d >> log'%(d,i,i,args.t))
        with open('msa/blast.out-%d'%i) as fin:
            record = fin.readline().strip().split()
            for seq_record in SeqIO.parse('%s/assembly-%d.fa'%(d,i), 'fasta'):
                if int(record[8]) < int(record[9]):
                    fout.write(str(seq_record.seq)+'\n')
                else:
                    fout.write(str(seq_record.seq.reverse_complement())+'\n')

print('Running MAFFT ...\n')
os.system('cd msa; mafft --thread %d msaInput > msaOutput 2>/dev/null'%args.t)

print('Generating the final assembly and computing the entropy ...\n')
multiSeq=[]
entropy=[]
consensus=[]
for seq_record in SeqIO.parse('msa/msaOutput', 'fasta'):
    multiSeq.append(list(str(seq_record.seq)))
for i in range(len(multiSeq[0])):
    candidate=[]
    for seq in multiSeq:
        candidate.append(seq[i])
    candidate=collections.Counter(candidate)
    e=0
    for f in candidate.values():
        p=f/sum(candidate.values())
        if p>0:
            e+=-p*np.log(p)
    entropy.append(e)
    base=max(candidate,key=lambda k:candidate[k]).upper()
    if base=='-':
        base=''
    consensus.append(base)

start=-1
end=len(consensus)
for i in range(len(entropy)):
    if entropy[i]>0.5:
        if i<500:
            start=i
        elif i>len(consensus)-500:
            end=i
            break

consensus=consensus[start+1:end]
entropy=entropy[start+1:end]

start=0
while consensus[start]=='':
    start+=1
end=len(consensus)-1
while consensus[end]=='':
    end-=1
consensus=consensus[start:end+1]
entropy=entropy[start:end+1]

consensus=''.join(consensus)
maxEntropy=max(entropy)

with open('finalAssembly.fa','w') as fout:
    fout.write('>contig1\n')
    fout.write(consensus+'\n')

# with open('entropy.perBase', 'w') as fout:
#     fout.write('>entropy1\n')
#     for k in range(len(entropy)):
#         if (k+1)%100 == 0:
#             fout.write(str(round(entropy[k],4))+'\n')
#         else:
#             fout.write(str(round(entropy[k],4))+' ')
#     fout.write('\n')
    
if maxEntropy>0.5:
    print('The assembly is not stable and its entropy is %.4f.\n\nYou may change the sampling coverage or tune the parameters.'%maxEntropy)
else:
    print('\nFinished! The final assembly is stable and its entropy is %.4f.\n'%maxEntropy)
