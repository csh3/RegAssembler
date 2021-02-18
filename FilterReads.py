# Copyright © 2021, Shenghao Cao & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-r1', required=True, help='Fastq file with forward paired reads')
parser.add_argument('-r2', required=True, help='Fastq file with reverse paired reads')
parser.add_argument('-l', required=True, type=int, help='Output reads no shorter than this value')
parser.add_argument('-q', default=0, type=int, help='Output reads with all the Phred scores no less than this value')
parser.add_argument('-o1', default='filteredReads1.fq', help='Output file with forward filtered reads')
parser.add_argument('-o2', default='filteredReads2.fq', help='Output file with reverse filtered reads')

args = parser.parse_args()

record_iterator_1 = SeqIO.parse(args.r1, "fastq")
record_iterator_2 = SeqIO.parse(args.r2, "fastq")
length = args.l
quality = args.q

filtered_reads_1 = []
filtered_reads_2 = []

while True:
    try:
        rec_1 = next(record_iterator_1)
        rec_2 = next(record_iterator_2)
        if len(rec_1.seq)>=length and len(rec_2.seq)>=length and min(rec_1.letter_annotations["phred_quality"]) >= quality and min(rec_2.letter_annotations["phred_quality"]) >= quality:
            filtered_reads_1.append(rec_1)
            filtered_reads_2.append(rec_2)
    except StopIteration:
        break

count1=SeqIO.write(filtered_reads_1, args.o1, "fastq")
count2=SeqIO.write(filtered_reads_2, args.o2, "fastq")

