# Copyright © 2021, Shenghao Cao & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

# from Bio import SeqIO
import sys
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-r1', required=True, help='Fastq file with forward paired reads')
parser.add_argument('-r2', required=True, help='Fastq file with reverse paired reads')
parser.add_argument('-l', required=True, type=int, help='Output reads no shorter than this value')
# parser.add_argument('-q', default=0, type=int, help='Output reads with all the Phred scores no less than this value')
parser.add_argument('-o1', default='filteredReads1.fq', help='Output file with forward filtered reads')
parser.add_argument('-o2', default='filteredReads2.fq', help='Output file with reverse filtered reads')

args = parser.parse_args()
length = args.l

readFile_1 = open(args.r1,'r')
readFile_2 = open(args.r2,'r')
total_reads_1 = readFile_1.readlines()
total_reads_2 = readFile_2.readlines()
readFile_1.close()
readFile_2.close()

with open("temp-1.fq",'w') as fout1:
    with open("temp-2.fq",'w') as fout2:
        for k in range(int(len(total_reads_1)/4)):
            if len(total_reads_1[4*k+1])-1 >= length and len(total_reads_2[4*k+1])-1 >= length:
                fout1.write(total_reads_1[4*k])
                fout1.write(total_reads_1[4*k+1])
                fout1.write(total_reads_1[4*k+2])
                fout1.write(total_reads_1[4*k+3])
                fout2.write(total_reads_2[4*k])
                fout2.write(total_reads_2[4*k+1])
                fout2.write(total_reads_2[4*k+2])
                fout2.write(total_reads_2[4*k+3])

with open('inputFile','w') as fout:
	fout.write('temp-1.fq\n')
	fout.write('temp-2.fq\n')

os.system('fastuniq -i inputFile -t q -o %s -p %s'%(args.o1,args.o2))
os.system('rm inputFile temp-1.fq temp-2.fq')
# record_iterator_1 = SeqIO.parse(args.r1, "fastq")
# record_iterator_2 = SeqIO.parse(args.r2, "fastq")
# quality = args.q

# filtered_reads_1 = []
# filtered_reads_2 = []

# while True:
#     try:
#         rec_1 = next(record_iterator_1)
#         rec_2 = next(record_iterator_2)
#         if len(rec_1.seq)>=length and len(rec_2.seq)>=length and min(rec_1.letter_annotations["phred_quality"]) >= quality and min(rec_2.letter_annotations["phred_quality"]) >= quality:
#             filtered_reads_1.append(rec_1)
#             filtered_reads_2.append(rec_2)
#     except StopIteration:
#         break

# count1=SeqIO.write(filtered_reads_1, args.o1, "fastq")
# count2=SeqIO.write(filtered_reads_2, args.o2, "fastq")

