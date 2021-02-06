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
3. `BWA (version 0.7.17-r1188)`
4. `MAFFT (v7.467)`
5. `SPAdes (v3.14.0)` This is optional if the user would like to try a SPAdes version of the re-sampling scheme.

The latter three tools can be installed using Bioconda (https://bioconda.github.io/). 

**Please ensure the above in your environment variable `PATH`.**

## 4. Commands
The main program is `Re-sampling.py`, and we recommand the following parameters.

For **Miseq** datasets:

```
python Re-sampling.py -thr 10 -ho 5 -al 30 -t THREADS -s SAMPLES
```

For datasets sequenced from other illumina platforms (such as **NextSeq 550** and **NovaSeq 6000**):

```
python Re-sampling.py -thr 3 -ho 2 -al 20 -t THREADS -s SAMPLES
```

For **simulated** datasets:

```
python Re-sampling.py -nchi -thr 3 -ho 2 -al 30 -t THREADS -s SAMPLES
```

`-t THREADS` specifies the number of threads you would like to use. 

`-s SAMPLES` specifies the number of samplings in the re-sampling scheme. The default setting is `-s 10`, but you can specify `-s 1` to disuse the re-sampling scheme.

## 5. Current version

The version of the current release is v1.0.


## 6. Contact

Please contact `cao.shenghao@foxmail.com` for any questions.


## 7. License

**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License**

For details, please read `RegAssembler/License`.
