Copyright © 2021, Shenghao Cao & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

# RegAssembler

## 1. Introduction
RegAssembler is a genome assembler employing the robust regression and re-sampling techniques.

The current version is specially designed for high-confidence reconstruction of SARS-CoV-2 genomes using Illumina sequencing reads.

## 2. Installation
You can download the software package by the command:

```
git clone https://github.com/csh3/RegAssembler.git
```

Then run the pipeline in the directory `RegAssembler`.

## 3. Dependencies
The current version requires:

1. `Python3` with the following modules: 
`os`, `sys`, `re`, `argparse`, `biopython`, `numpy`, `math`, `random`, `networkx`, `scipy`, `statsmodels`, `copy`, `collections`, `time`, `multiprocessing`

2. `BLAST (version 2.9.0+)`
3. `BWA (version 0.7.17)`
4. `MAFFT (v7.467)`
5. `fastuniq (version 1.1-1)`
6. `SPAdes (v3.14.0)` This is optional if the user would like to try the re-sampling version of SPAdes.

The latter five tools can be installed through the [Bioconda](https://bioconda.github.io/) channel. 

**Please ensure the above in your environment variable `PATH`.**

## 4. Usage
We recommend following parameter settings for different types of sequencing data. The main program is `Re-sampling.py` and `FilterReads.py` is used to remove short and duplicate reads in advance. 

For **MiSeq** datasets:

```
python FilterReads.py -r1 readFile_1 -r2 readFile_2 -l 300 -o1 filteredReads1.fq -o2 filteredReads2.fq
python Re-sampling.py -r1 filteredReads1.fq -r2 filteredReads2.fq -n1 10000 -n2 20000 -thr 10 -ho 5 -al 30 -t THREADS -s SAMPLES
```

For **NextSeq 550** datasets:

```
python FilterReads.py -r1 readFile_1 -r2 readFile_2 -l 150 -o1 filteredReads1.fq -o2 filteredReads2.fq
python Re-sampling.py -r1 filteredReads1.fq -r2 filteredReads2.fq -n1 15000 -n2 60000 -thr 3 -ho 2 -al 20 -t THREADS -s SAMPLES
```

For **simulated** datasets:

```
python Re-sampling.py -r1 readFile_1 -r2 readFile_2 -n1 10000 -n2 10000 -thr 3 -ho 2 -al 30 -nchi -t THREADS -s SAMPLES
```

`-t THREADS` specifies the number of CPU threads you would like to use. 

`-s SAMPLES` specifies the number of samplings. The default setting is `-s 10`, but you can specify `-s 1` to disuse the re-sampling scheme.

## 5. Example
We offer an example that runs RegAssembler on a MiSeq dataset. First download the dataset from NCBI's Sequence Read Archive using SRA Toolkit.

```
prefetch SRR12089766
fastq-dump --split-e SRR12089766
```

Then filter out short reads.

```
python FilterReads.py -r1 SRR12089766_1.fastq -r2 SRR12089766_2.fastq -l 300 -o1 filteredReads1.fq -o2 filteredReads2.fq
```
Run RegAssembler with 1 sampling.

```
python Re-sampling.py -r1 filteredReads1.fq -r2 filteredReads2.fq -n1 10000 -n2 20000 -thr 10 -ho 5 -al 30 -s 1 -t THREADS
```
Run RegAssembler with 10 samplings.

```
python Re-sampling.py -r1 filteredReads1.fq -r2 filteredReads2.fq -n1 10000 -n2 20000 -thr 10 -ho 5 -al 30 -s 10 -t THREADS
```
`-t THREADS` specifies the number of CPU threads you would like to use. 

## 6. Current version

The version of the current release is v1.0.


## 7. Contact

Please contact <cao.shenghao@foxmail.com> for any questions.

## 8. License

**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License**

For details, please read `RegAssembler/License`.
