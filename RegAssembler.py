# Copyright © 2021, Shenghao Cao & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

import os
import numpy as np
import networkx as nx
import scipy.sparse as sp
import copy
import collections
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align
import argparse
import time
from multiprocessing import Pool,Manager


def rename_input(inp1, inp2):

    reads=[]
    number=-1
    with open('trainReads.fasta','w') as fout:
        for seq_record in SeqIO.parse(inp1, 'fastq'):
            number+=1
            fout.write('>'+str(number)+'\n')
            fout.write(str(seq_record.seq)+'\n')
            reads.append(str(seq_record.seq))

        for seq_record in SeqIO.parse(inp2, 'fastq'):
            number+=1
            fout.write('>'+str(number)+'\n')
            fout.write(str(seq_record.seq)+'\n')
            reads.append(str(seq_record.seq))

    readsNum=len(reads)

    return(readsNum, reads)


def detectChimeras(blast, chimeras_dic, thread, chi=15, threshold=2):
    chimeras_list = []
    read_tmp = ''
    highScore_num = np.infty
    
    with open(blast) as fin:
        for line in fin:
            record = line.strip().split()
            if record[0] == read_tmp:
                if record[0] != record[1] and int(record[3])>len(reads[int(read_tmp)])-chi and not (int(record[6])==int(record[9]) and int(record[7])==int(record[8])):
                    highScore_num += 1
            else:
                if highScore_num < threshold:
                    chimeras_list.append(int(read_tmp))
                read_tmp = record[0]
                highScore_num = 0
                if record[0] != record[1] and int(record[3])>len(reads[int(read_tmp)])-chi and not (int(record[6])==int(record[9]) and int(record[7])==int(record[8])):
                    highScore_num += 1
        if highScore_num < threshold:
            chimeras_list.append(int(read_tmp))
        
    chimeras_dic[thread] = chimeras_list


def removeChimeras(reads):
    chimeras_list = []
    for thread in chimeras_dic:
        chimeras_list += chimeras_dic[thread]
    index = sorted(list(set(range(len(reads)))-set(chimeras_list)))
    reads = [reads[k] for k in index]
    readsNum=len(reads)
    for i in range(readsNum):
        reads.append(str(Seq(reads[i]).reverse_complement()))

    return(reads, readsNum, np.array(chimeras_list))


def setupRegressionModel(blast, model_dic, thread, hangingOut=2, alignmentLen=20):
    
    n=0
    
    row=0
    indptr=[0]
    indices=[]
    data=[]
    response=[]


    with open(blast) as fin:
        for line in fin:

            record = line.strip().split()
            i = int(record[0])
            j = int(record[1])
            len1 = len(reads[i])
            len2 = len(reads[j])
            if i<j and i not in chimeras_list and j not in chimeras_list:
                if int(record[8])<int(record[9]) and min(int(record[6])-1, int(record[8])-1)<=hangingOut and min(len1-int(record[7]), len2-int(record[9]))<=hangingOut and int(record[3]) >= alignmentLen:
                    row+=2
                    
                    i -= len(np.where(chimeras_list<i)[0])
                    j -= len(np.where(chimeras_list<j)[0])

                    y=len2-int(record[9])-(len1-int(record[7]))
                    response.append([y])
                    indices.append(j)
                    data.append(1)
                    indices.append(i)
                    data.append(-1)
                    indptr.append(indptr[-1]+2)

                    y=(int(record[8])-1)-(int(record[6])-1)
                    response.append([y])
                    indices.append(i+readsNum)
                    data.append(-1)
                    indices.append(j+readsNum)
                    data.append(1)
                    indptr.append(indptr[-1]+2)

                elif int(record[8])>int(record[9]) and min(int(record[6])-1, len2-int(record[8]))<=hangingOut and min(len1-int(record[7]), int(record[9])-1)<=hangingOut and int(record[3]) >= alignmentLen:
                    row+=2
                    
                    i -= len(np.where(chimeras_list<i)[0])
                    j -= len(np.where(chimeras_list<j)[0])

                    y=int(record[9])-1-(len1-int(record[7]))
                    response.append([y])
                    indices.append(j+readsNum)
                    data.append(1)
                    indices.append(i)
                    data.append(-1)
                    indptr.append(indptr[-1]+2)

                    y=int(record[6])-1-(len2-int(record[8]))
                    response.append([y])
                    indices.append(i+readsNum)
                    data.append(1)           
                    indices.append(j)
                    data.append(-1)
                    indptr.append(indptr[-1]+2)
                               

    design=sp.csr_matrix((data, indices, indptr), shape=(row, len(reads)))
    if response==[]:
        response=np.matrix([[]]).T
    else:
        response=np.matrix(response)

    model_dic[thread]=[design, response]


def deleteRowsCsr(mat, indices):
    indices = list(indices)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask]


def IRLS(X,Y,THR=3,ridge=False):
    t=X.T
    A=t.dot(X)
    y=sp.csr_matrix(Y)
    b=t.dot(y).todense()
    if ridge:
        estimate=sp.linalg.spsolve(A+0.001*sp.identity(A.shape[0]),b)
    else:
        estimate, exitCode = sp.linalg.lgmres(A, b, atol=1e-05)
    residual=abs((X.dot(sp.csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    max_residual=max(residual)
    print("max_residual:", round(max_residual,2))

    while max_residual>THR:
        if max_residual>1000:
            threshold=1000
        elif max_residual>100:
            threshold=100
        elif max_residual>10:
            threshold=10
        elif max_residual>5:
            threshold=5
        else:
            threshold=max_residual-1


        residual_odd = residual[1::2]
        residual_even = residual[0::2]
        residual_max=np.zeros(len(residual_odd))
        index = residual_odd>residual_even
        residual_max[index] = residual_odd[index]
        residual_max[~index] = residual_even[~index]
        residual=np.zeros(len(residual))
        residual[0::2]=residual_max
        residual[1::2]=residual_max

        index=np.where(residual>threshold)[0]
        residual=np.delete(residual,index)
        reweight=np.square(1-np.square(residual/threshold))
        # reweight=np.ones(len(residual))
        reweight=sp.diags(reweight)
        X=deleteRowsCsr(X,index)
        Y=np.delete(Y,index,0)
        t=X.T
        A=t.dot(reweight).dot(X)
        y=sp.csr_matrix(Y)
        b=t.dot(reweight).dot(y).todense()
        if ridge:
            estimate=sp.linalg.spsolve(A+0.001*sp.identity(A.shape[0]),b)
        else:
            estimate, exitCode = sp.linalg.lgmres(A, b, atol=1e-05)
        residual=abs((X.dot(sp.csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        max_residual=max(residual)
        print("max_residual:", round(max_residual,2))

    G = nx.Graph()
    for i in range(int(len(X.indices)/2)):
        n1=X.indices[2*i]
        n2=X.indices[2*i+1]
        G.add_edge(n1,n2)

    reads_list=[]
    estimates_list=[]
    index_list=[]
    for c in nx.connected_components(G):
        sub_index=list(G.subgraph(c).nodes)
        sub_estimates=[estimate[i] for i in sub_index]
        sub_reads=[reads[i] for i in sub_index]
        order=np.argsort(sub_estimates)
        sub_index=[sub_index[i] for i in order]
        sub_reads=[sub_reads[i] for i in order]
        sub_estimates=[sub_estimates[i] for i in order]
        if (sub_index[0]+readsNum not in sub_index) and (sub_index[0]-readsNum not in sub_index):
            loc=len(index_list)
            for k in range(len(index_list)):
                if (sub_index[0]+readsNum in index_list[k]) or (sub_index[0]-readsNum in index_list[k]):
                    loc=k+1
            index_list.insert(loc, sub_index)
            reads_list.insert(loc, sub_reads)
            estimates_list.insert(loc, sub_estimates)

    return(estimates_list, reads_list, index_list)


def linkContigs(estimates_list, reads_list, index_list, X):

    fragments=list(range(int(len(estimates_list)/2)))
    Estimates_list=[]
    Reads_list=[]
    A=X.T.dot(X)

    while len(fragments)>0:
        fragInitial=fragments.pop(0)*2
        Estimates_list.append([estimates_list[fragInitial]])
        Reads_list.append([reads_list[fragInitial]])
        for i in range(2):
            frag=fragInitial
            while True:
                read=reads_list[frag]
                index=index_list[frag]
                estimate=estimates_list[frag]

                if i==0:
                    ind_list=range(len(estimate))[-25:]
                else:
                    ind_list=range(len(estimate))[:25]
                overlapLink=[]
                bridge={}
                for ind in ind_list:
                    a=A[:,index[ind]]
                    for j in [2*k for k in fragments]+[2*k+1 for k in fragments]:
                        for column in a.indices:
                            if (i==0 and column in index_list[j][:25]) or (i==1 and column in index_list[j][-25:]):
                                overlapLink.append(j)
                                bridge[j]=(ind,index_list[j].index(column))
                            
                overlapCount=collections.Counter(overlapLink)
                success=False
                for successor in sorted(list(overlapCount.keys()), key=lambda k:overlapCount[k], reverse=True):
                    estimate_s=estimates_list[successor]
                    read_s=reads_list[successor]
                    alignments = aligner.align(read[bridge[successor][0]], read_s[bridge[successor][1]])
                    alignment = alignments[0]
                    if alignment.score>=20 and overlapCount[successor]>=10:
                        # print(alignment.score, (frag,successor))
                        translocation=estimate[bridge[successor][0]]-(len(read[bridge[successor][0]])-alignment.aligned[0][-1][1])+(len(read_s[bridge[successor][1]])-alignment.aligned[1][-1][1])-estimate_s[bridge[successor][1]]
                        estimates_list[successor]=[k+translocation for k in estimate_s]
                        Estimates_list[-1].append(estimates_list[successor])
                        Reads_list[-1].append(read_s)
                        frag=successor
                        fragments.remove(int(successor/2))
                        success=True
                        break
                if not success:
                    break

    for k in range(len(Estimates_list)):
        Estimates_list[k]=[j for i in Estimates_list[k] for j in i]
        Reads_list[k]=[j for i in Reads_list[k] for j in i]

    return(Estimates_list, Reads_list)


def generateReference(reads, estimates):

    coordinates=[estimates[i]-len(reads[i]) for i in range(len(estimates))]
    length = max(estimates)-min(coordinates)

    reference_element=[dict() for i in range(length)]
    for i in range(len(reads)):
        read = reads[i]
        distance=estimates[i]-len(read)-min(coordinates)
        for j in range(len(read)):
            if read[j] in reference_element[j+distance]:
                reference_element[j+distance][read[j]]+=1
            else:
                reference_element[j+distance][read[j]]=1
    for k in range(length):
        if sum(reference_element[k].values())<1:
            reference_element[k]=''
        else:
            reference_element[k] = max(reference_element[k], key=reference_element[k].get)
    reference=''.join(reference_element)

    return(reference)


def remapToReference(reference, reads, estimates, minCoordinate):

    consensus=[]
    for base in reference:
        consensus.append({base:0, 'insertion':{'':0}})

    for k in range(len(reads)):
        query=reads[k]
        location=estimates[k]-minCoordinate
        startPoint=max(location-len(query)-20,0)
        ref=reference[startPoint:location+20]
        alignments = aligner.align(query, ref)
        alignment = alignments[0]
        if alignment.score < 20:
            continue
        cigar_query = alignment.aligned[0]
        cigar_ref = alignment.aligned[1]
        for i in range(len(cigar_ref)):
            L=0
            R=cigar_ref[i][1]-cigar_ref[i][0]
            if i==0:
                L=3
            if i==len(cigar_ref)-1:
                R=R-3
            for j in range(L,R):
                base=query[cigar_query[i][0]+j]
                if base in consensus[startPoint+cigar_ref[i][0]+j]:
                    consensus[startPoint+cigar_ref[i][0]+j][base]+=1
                else:
                    consensus[startPoint+cigar_ref[i][0]+j][base]=1
                if j>0:
                    consensus[startPoint+cigar_ref[i][0]+j]['insertion']['']+=1

        for i in range(len(cigar_ref)-1):
            left=cigar_ref[i][1]
            right=cigar_ref[i+1][0]
            if right-left > 0:
                for j in range(startPoint+left, startPoint+right):
                    if '' in consensus[j]:
                        consensus[j]['']+=1
                    else:
                        consensus[j]['']=1
                    consensus[j]['insertion']['']+=1
                consensus[startPoint+right]['insertion']['']+=1
            else:
                insertion=query[cigar_query[i][1]:cigar_query[i+1][0]]
                if insertion in consensus[startPoint+left]['insertion']:
                    consensus[startPoint+left]['insertion'][insertion]+=1
                else:
                    consensus[startPoint+left]['insertion'][insertion]=1
        # if (k+1)%1000==0:
        #     print('thread '+str(thread)+' mapped '+str(k+1)+' reads ')

    for k in range(len(consensus)):
        insertion_dic=consensus[k].pop('insertion')
        insertion=max(insertion_dic, key=insertion_dic.get)
        base=max(consensus[k], key=consensus[k].get)
        consensus[k]=insertion+base

    consensus=('').join(consensus)
    return(consensus)


def joinContigs(consensus_list, minScore=20, trimLen=500):
    contigs = sorted(consensus_list, key=lambda k:len(k), reverse=True)
    contig_list = []
    index = list(range(len(contigs)))
    mergeNum=0
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
                mergeNum+=1
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
                mergeNum+=1
                ref = maxQuery[:-trimLen]+maxQuery[-trimLen:][:maxAlignment.aligned[0][-1][-1]]+ref[:trimLen][maxAlignment.aligned[1][-1][-1]:]+ref[trimLen:]
                index.remove(maxIndex)
            else:
                break
        if len(ref) > args.cl:
            contig_list.append(ref)
    print("\nJoin %d consensus"%mergeNum)
    return(contig_list)


def writeResults(contig_list,outFile):
    with open(outFile,'w') as fout:
        n=0
        for contig in contig_list:
            n+=1
            fout.write('>contig'+str(n)+'\n')
            line = int(len(contig)/100)
            for k in range(line):
                fout.write(contig[100*k:100*(k+1)]+'\n')
            if len(contig) > line*100:
                fout.write(contig[line*100:]+'\n')

    table = "{0:<10}\t{1:<10}"
    print(table.format('contigName', 'length'))
    contigNum=0
    for contig in contig_list:
        contigNum += 1
        print(table.format('contig%d'%contigNum, len(contig)))
    print(table.format('\ntotal', sum([len(k) for k in contig_list])))


if __name__ == "__main__":

    descript="This program assembles genomes by robust regression.\n"
    parser = argparse.ArgumentParser(description=descript)
    parser.add_argument('-r1', required=True, help='Fastq file with forward paired reads')
    parser.add_argument('-r2', required=True, help='Fastq file with reverse paired reads')
    parser.add_argument('-o', default='draft.fa', help='Output file')
    parser.add_argument('-t', type=int, default=1, help='Number of threads for parallelism')
    parser.add_argument('-thr', type=int, default=3, help='Residual threshold for IRLS algorithm to halt')
    parser.add_argument('-ho', type=int, default=2, help='Admissible hanging-out length for pairwise overlaps')
    parser.add_argument('-ip', type=int, default=98, help='Minimum identity percentage for a successful overlap')
    parser.add_argument('-al', type=int, default=20, help='Minimum alignment length for a successful overlap')
    parser.add_argument('-nchi', action="store_true", help='No chimeric reads need to be detected and removed')
    parser.add_argument('-cl', type=int, default=500, help='Shorter contigs will be removed')
    # parser.add_argument("-ridge", action="store_true", help="Use robust ridge regression")
    args = parser.parse_args()

    print("\n--------------------------------------\n")
    print("Parameters: RegAssembler -r1 %s -r2 %s -o %s -t %d -thr %d -ho %d -ip %d -al %d -cl %d -nchi %s \n"%(args.r1,args.r2,args.o,args.t,args.thr,args.ho,args.ip,args.al,args.cl,args.nchi))

    startTime = time.time()

    aligner = Align.PairwiseAligner(mode = 'local', match_score=1, mismatch_score=-2, open_gap_score=-5, extend_gap_score=-1)

    print("\n--------------------------------------\nDetect overlaps among reads by Blastn.\n")

    readsNum, reads= rename_input(args.r1, args.r2)
    recordNum=int(np.ceil(readsNum/args.t))
    print('\nMaking  Blastn database...\n')
    os.system("makeblastdb -in trainReads.fasta -dbtype nucl -out db > /dev/null 2>&1")
    if os.path.exists("blastnResults"):
        os.system("rm -rf blastnResults")
    os.system("mkdir blastnResults")
    print('\nRunning Blastn...')
    os.system("cat trainReads.fasta | parallel --recstart '>' -N %d --pipe blastn -evalue 0.01 -gapopen 5 -gapextend 2 -penalty -2 -reward 1 -word_size 10 -perc_identity %d -outfmt 6 -db db -query - -out 'blastnResults/blastnResult-{#}'"%(recordNum, args.ip))

    m=Manager()

    if not args.nchi:
        print("\n--------------------------------------\nDetect and remove potential chimeras.\n")
        pool1 = Pool(args.t)
        chimeras_dic = m.dict()
        for thread in range(args.t):
            thread+=1
            pool1.apply_async(func=detectChimeras, args=('blastnResults/blastnResult-'+str(thread), chimeras_dic, thread))
        pool1.close()
        pool1.join()
        reads, readsNum, chimeras_list = removeChimeras(reads)
    else:
        for i in range(readsNum):
            reads.append(str(Seq(reads[i]).reverse_complement()))
        chimeras_list = np.array([])

    print("\n--------------------------------------\nSolve coordinates of reads by regression.\n")

    print("\nScanning records of overlapping reads...\n")
    pool1 = Pool(args.t)
    model_dic=m.dict()
    for thread in range(args.t):
        thread+=1
        pool1.apply_async(func=setupRegressionModel,args=('blastnResults/blastnResult-'+str(thread), model_dic, thread, args.ho, args.al))
    pool1.close()
    pool1.join()

    design_list=[]
    response_list=[]
    for thread in range(args.t):
        thread+=1
        design_list.append(model_dic[thread][0])
        response_list.append(model_dic[thread][1])

    design=sp.vstack(design_list)
    response=np.vstack(response_list)

    print("\nRunning iteratively reweighted least squares algorithm...\n")
    estimates_list, reads_list, index_list= IRLS(copy.deepcopy(design),copy.deepcopy(response),args.thr)

    print("\nLinking separate blocks...")
    Estimates_list, Reads_list = linkContigs(estimates_list, reads_list, index_list, copy.deepcopy(design))

    print("\n--------------------------------------\nGenerate consensus by coordinates of reads.\n")
    consensus_list=[]
    contigNum=0
    for i in range(len(Estimates_list)):
        read=Reads_list[i]
        estimate=list(map(int,np.round(Estimates_list[i])))
        reference=generateReference(read, estimate)
        if len(reference) > 500:
            contigNum+=1
            minCoordinate=min([estimate[i]-len(read[i]) for i in range(len(read))])

            print('Generating consensus%d ...'%contigNum)
            consensus=remapToReference(reference, read, estimate, minCoordinate)

            consensus_list.append(consensus)

    contig_list = joinContigs(consensus_list)

    print("\n--------------------------------------\nWrite results.\n")

    writeResults(contig_list,args.o)

    endTime = time.time()
    executionTime = endTime-startTime
    hour = int(executionTime/3600)
    minute = int((executionTime-hour*3600)/60)
    second = round(executionTime-hour*3600-minute*60,2)

    print("\n--------------------------------------\nFinish draft assembly.\n")
    print("\nExecution time: "+(str(hour)+"h")*bool(hour)+(str(minute)+"m")*bool(minute)+(str(second)+"s")*bool(second)+"\n")