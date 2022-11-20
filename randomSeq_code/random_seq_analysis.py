from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce

def get_fasta_some_seqs(fastaLis, nSeq, sel_idx, fnm):
    # read in FASTA files and rev. complement for positive and negatives
    fout = open(fnm, 'w')
    for f_idx, f in enumerate(fastaLis):
        base_ = f_idx*nSeq
        print(base_)
        for idx, seq in enumerate(SeqIO.parse(f, "fasta")):
            if base_+ idx in sel_idx:
                print(seq.seq, file = fout)
            elif idx == nSeq:
                break
    fout.close()
    return 0

fnm_random = ['/home/mleone2/ml_class_project/results/random_seqs_test' + str(i) + '_predictions_10000result.txt' for i in range(1, 6)]

arr_random = np.array([])
for fnm in fnm_random:
    with open(fnm, 'r') as inp:
        for line in inp:
            arr_random = np.append(arr_random, float(line.strip()))

idx_high = np.where(arr_random>0.98)[0]
idx_low = np.where(arr_random<0.00001)[0]

fasta_random = ['/home/mleone2/ml_class_project/results/random_seqs_test' + str(i) + '.fa' for i in range(1, 6)]
f_high = '/home/mleone2/ml_class_project/results/high_score.txt'
f_low = '/home/mleone2/ml_class_project/results/low_score.txt'
high_seq = get_fasta_some_seqs(fasta_random, 10000, idx_high, f_high)
low_seq = get_fasta_some_seqs(fasta_random, 10000, idx_low, f_low)

# K = 5
# def kmer_cnt(seq, K, d):
#     for i in range(len(seq) - K + 1):
#         d[seq[i:i+K]] += 1
#     return d
