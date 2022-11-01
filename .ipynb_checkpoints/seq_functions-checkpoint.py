# from tensorflow.keras.models import Sequential
# from tensorflow.keras.layers import Dense, Dropout, Activation, Embedding, Conv2D, MaxPooling2D, Flatten
# from tensorflow.keras.optimizers import SGD, Adam
# from tensorflow.keras.callbacks import EarlyStopping
# from tensorflow.keras.regularizers import l2
from tensorflow.keras.models import load_model
#from clr_callback import *
from sklearn import metrics
from Bio import SeqIO

import tensorflow as tf
# import tensorflow.keras.backend as K
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import math
import os

import time

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
# make sure running on GPU #
tf.config.experimental.list_physical_devices('GPU')
# tf.get_logger().setLevel('WARNING')

def onehot_seq(seq):
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    to_return = np.zeros((len(seq),4), dtype='int8')
    for idx,letter in enumerate(seq):
        if letter not in ['N','n']:
            to_return[idx,letter_to_index[letter]] = 1
    return to_return

def inverse_onehot_seq(seq):
    ### input: (length,4)
    ### output: list, being ['A','C','G','T'] of bases
    
    letter_to_index =  {'A':0,
                        'C':1,
                        'G':2,
                        'T':3}
    index_to_letter = {v: k for k, v in letter_to_index.items()}
    
    seq_as_letters = []
    for base_dex in range(seq.shape[0]):
        letter = index_to_letter[   np.nonzero(seq[base_dex])[0][0]   ]
        seq_as_letters.append(letter)
        
    return seq_as_letters

def get_fasta_seqs(fasta):
    # read in FASTA files and rev. complement for positive and negatives
    seqs = np.array([onehot_seq(seq) for seq in SeqIO.parse(fasta, "fasta") ])
    
    return seqs

def get_fasta_some_seqs(fasta, nSeq, lenSeq):
    # read in FASTA files and rev. complement for positive and negatives
    seqs = np.zeros((nSeq, lenSeq, 4), dtype='int8')
    for idx, seq in enumerate(SeqIO.parse(fasta, "fasta")):
        if idx == nSeq:
            break
        else:
            seqs[idx] = onehot_seq(seq)
    return seqs
# def encode_sequence(fasta_pos, fasta_neg, shuffleOff = True):
#     # read in FASTA files and rev. complement for positive and negatives
#     x_pos = np.array([onehot_seq(seq) for seq in SeqIO.parse(fasta_pos, "fasta") ] +
#     [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(fasta_pos, "fasta") ]) 
#     x_neg = np.array([onehot_seq(seq) for seq in SeqIO.parse(fasta_neg, "fasta") ] +
#     [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(fasta_neg, "fasta") ])
#     # concatenate positives and negatives
#     print(f'There {x_pos.shape[0]} positives and {x_neg.shape[0]} negatives.')
#     x = np.expand_dims(np.concatenate((x_pos, x_neg)), axis=3) 
#     y = np.concatenate((np.ones(len(x_pos)),np.zeros(len(x_neg))))
#     # need to shuffle order of training set for validation splitting last
#     if not shuffleOff:
#         indices = np.arange(y.shape[0])
#         np.random.shuffle(indices)
#         x = x[indices,:]
#         y = y[indices]
#     #
#     return x, y

def quick_predict_sequences(model, seqs):
    
    y_pred_score = model.predict(seqs)
    return y_pred_score

def main():
    print('this is a test demonstrating the functions and how to save predictions')
    model_name = '/home/mleone2/ml_class_project/pretrained_models/Mo2015_EXCpos_Ctx_fold1/Mo2015_EXCpos_Ctx_fold1_OCP_NB1000_NE23_BR0.01_MR0.1_BM0.85_MM0.99_DO0.1.h5'
    parser = argparse.ArgumentParser(description='predictions')
    parser.add_argument("--fasta_name")
    parser.add_argument("--prediction_file")
    args = parser.parse_args()
    
    # path = '/home/mleone2/ml_class_project/results/'
    # fasta_name = 'random_seqs_test'
    # full_fasta = path + fasta_name + '.fa' 
    full_fasta = args.fasta_name
    prediction_file = args.prediction_file

    t0 = time.time()
    print('load 100k sequences')
    # seqs = get_fasta_seqs(full_fasta)
    some_seqs = get_fasta_some_seqs(full_fasta, 1000, 501)

    t1 = time.time()
    print(str(t1-t0) + ' seconds')
    
    # some_seqs = seqs[0:1000,:,:]
    
    print('load model then predict for 100 sequences')
    print('loading model')
    model = load_model(model_name, compile=False)
    print('predict for 100 sequences')
    predictions = quick_predict_sequences(model_name, some_seqs)
    t2 = time.time()
    print(str(t2-t1) + ' seconds')
    
    print('save predictions')
    # prediction_file = path + fasta_name + '_predictions' + '.txt'
    np.savetxt(prediction_file, predictions)
    print('done')

if __name__ == "__main__":
    main()
    
    