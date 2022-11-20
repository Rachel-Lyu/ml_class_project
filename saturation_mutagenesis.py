#####
import seq_functions
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
import matplotlib.pyplot as plt

def saturation_mutate(seqs):
    
    # repeat sequences 4 times, for every base per position
    og_string_len = seqs.shape[1]
    # need 4 * og_string_len versions of each sequences
    seqs_repeated = np.repeat(seqs, 4*og_string_len, axis=0)
    
    ## og shift - shifting to build off
    
    og_shift = np.asarray([1+4*n*og_string_len for n in range(math.floor((len(seqs_repeated))/(4*og_string_len)))])
    
    #og_shift = np.asarray([4*n*og_string_len for n in range(math.floor((len(seqs_repeated))/(4*og_string_len)))])
    
    for string_id in range(og_string_len):
    
        one_shift = og_shift + 4*string_id
        two_shift = one_shift + 1
        three_shift = one_shift + 2
        
        #if string_id % 100 == 0:
        #    print("Onto base number " + str(string_id))
        
        seqs_repeated[one_shift,None,[string_id],:] = np.concatenate((np.expand_dims(seqs_repeated[one_shift,None,[string_id],-1],2), 
                                                   seqs_repeated[one_shift,None,[string_id],:-1]),2)
        
        seqs_repeated[two_shift,None,[string_id],:] = np.concatenate((seqs_repeated[two_shift,None,[string_id],-2:], 
                                                   seqs_repeated[two_shift,None,[string_id],:-2]),2)
        
        seqs_repeated[three_shift,None,[string_id],:] = np.concatenate((seqs_repeated[three_shift,None,[string_id],-3:], 
                                                   seqs_repeated[three_shift,None,[string_id],:-3]),2)

    return seqs_repeated

def saturation_mutagenesis(seqs,model, iterations, val_models = []):
    
    best_fitness_vals = []
    fitness = model.predict(seqs)
    best_fitness_vals.append(list(fitness))
    
    fitness_vals = []
    
    val_fitnesses_by_iter = []
    
    # predict initial validation fitnesses
    val_fitnesses = []
    if len(val_models) > 0:
        for val_model in val_models:
            val_fitness = val_model.predict(seqs)
            val_fitnesses.append(val_fitness)
    val_fitnesses = np.squeeze(  np.asarray(val_fitnesses)   )
    val_fitnesses_by_iter.append(val_fitnesses)
    
    ### temp
#     for seq in seqs:
#         print(seq_functions.inverse_onehot_seq(seq)[0:2]  )
    
    for step in range(iterations):
        print('')
        print('on step ' + str(step))
        seqs_mutated = saturation_mutate(seqs)
        
        ### temp
#         for seq_id in [0,1*501*4]:
#             for jj in [0,1,2,3,4,5]: #range(4*2):
#                 print(  seq_functions.inverse_onehot_seq(seqs_mutated[seq_id+jj,0:2,:])[0:2]  )
        
        
        t1 = time.time()
        fitness = model.predict(seqs_mutated)
        t2 = time.time()
        fitness_vals.append(fitness)
        
        print('avg fitness :' + str(np.mean(fitness)))
        
        seqs_to_keep = []
        
        best_fitness_this_round = []
        for seq_id in range(seqs.shape[0]):
            max_fitness = np.max(fitness[seq_id*4*seqs.shape[1]:(seq_id+1)*4*seqs.shape[1]])
            seq_to_keep = np.argmax(fitness[seq_id*4*seqs.shape[1]:(seq_id+1)*4*seqs.shape[1]])+seq_id*4*seqs.shape[1]
            
#             if seq_to_keep in seqs_to_keep:
#                 print(seq_to_keep)
#                 print(seqs_to_keep)
#                 print('found match - uh oh')
#                 exit()
            seqs_to_keep.append(seq_to_keep)
            
            best_fitness_this_round.append(max_fitness)
            
        assert len(seqs_to_keep) == seqs.shape[0]
        
        best_fitness_vals.append(best_fitness_this_round)
        new_seqs = seqs_mutated[seqs_to_keep,:,:]
        assert new_seqs.shape == seqs.shape
        seqs = new_seqs
        
        val_fitnesses = []
        if len(val_models) > 0:
            for val_model in val_models:
                val_fitness = val_model.predict(new_seqs)
                val_fitnesses.append(val_fitness)
        val_fitnesses = np.squeeze(  np.asarray(val_fitnesses)   )
        
        # each element of this list is a 4xnum_seq X 1
        val_fitnesses_by_iter.append(val_fitnesses)
        
        
    return seqs, np.transpose(np.asarray(best_fitness_vals)), np.asarray(fitness_vals), val_fitnesses_by_iter

def main():

    small_example_unonehot = False
    if small_example_unonehot:
        
        print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

        model_name = '/home/mleone2/ml_class_project/pretrained_models/Mo2015_EXCpos_Ctx_fold1/Mo2015_EXCpos_Ctx_fold1_OCP_NB1000_NE23_BR0.01_MR0.1_BM0.85_MM0.99_DO0.1.h5'
        path = '/home/mleone2/ml_class_project/results/'
        #fasta_name = 'random_seqs_test'
        #full_fasta = path + fasta_name + '.fa' 

        #some_seqs = seq_functions.get_fasta_some_seqs(full_fasta, 100, 501)
        
        my_seq_list = [['A']*501,['G']*501, ['C'] + ['T']*500]
        seqs = np.array([seq_functions.onehot_seq(seq) for seq in my_seq_list ])
        
        print('load model then predict for 3 sequences')
        print('loading model')
        model = load_model(model_name, compile=False)
        print('predict for 2 sequences')
        t1 = time.time()
        #predictions = seq_functions.quick_predict_sequences(model, some_seqs)
        t2 = time.time()
        print(str(t2-t1) + ' seconds')

        t3 = time.time()
        #some_seqs_mutated = saturation_mutate(some_seqs)
        iterations = 12
        seqs, best_fitness_vals, fitness_vals = saturation_mutagenesis(seqs,
                                                                       model,
                                                                       iterations)

        plt.figure()
        for ii in range(seqs.shape[0]):
            plt.plot(  range(iterations+1), best_fitness_vals[ii,:] )

        plt.savefig('fitness_small_test.png', bbox_inches = 'tight')

        np.save('all_fitnesses_smalltest.npy', fitness_vals, allow_pickle=True)
        np.save('best_fitnesses_smalltest.npy', best_fitness_vals, allow_pickle=True)

        t4 = time.time()
        #print(t4-t3)

        print('done')
        
    else:
        
        iterations = 15
        num_seqs = 200
        
        print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

        model_name = '/home/mleone2/ml_class_project/pretrained_models/Mo2015_EXCpos_Ctx_fold1/Mo2015_EXCpos_Ctx_fold1_OCP_NB1000_NE23_BR0.01_MR0.1_BM0.85_MM0.99_DO0.1.h5'
        path = '/home/mleone2/ml_class_project/results/'
        fasta_name = 'random_seqs_test'
        full_fasta = path + fasta_name + '.fa'
        
        val_model_names = [model_name.replace("fold1", "fold" + str(ii) ) for ii in [2,3,4,5] ]
        val_models = [load_model(name, compile=False) for name in val_model_names ]

        seqs_to_mutate = seq_functions.get_fasta_some_seqs(full_fasta, num_seqs, 501)
        np.save('../results/sat_mut_'+ str(num_seqs)+ 'seq_'+ str(iterations)+'iter_orig_seqs.npy', 
                seqs_to_mutate, allow_pickle=True)
        print('made orig seqs')
        
        #seqs_to_mutate = some_seqs#some_seqs[range(num_seqs),:,:]
        
        print('load model then predict for '+str(num_seqs)+' sequences')
        print('loading model')
        model = load_model(model_name, compile=False)
        print('predict for '+str(num_seqs)+' sequences')
        t1 = time.time()
        
        
        t2 = time.time()
        print(str(t2-t1) + ' seconds')
        t3 = time.time()
        
        seqs, best_fitness_vals, fitness_vals, val_fitnesses_by_iter = saturation_mutagenesis(seqs_to_mutate,model, iterations, val_models = val_models)

        plt.figure()
        for ii in range(seqs.shape[0]):
            plt.plot(  range(iterations+1), best_fitness_vals[ii,:] )

        np.save('../results/sat_mut_'+ str(num_seqs)+ 'seq_'+ str(iterations) + 'iter_all_fitnesses.npy', fitness_vals, allow_pickle=True)
        np.save('../results/sat_mut_'+ str(num_seqs)+ 'seq_'+ str(iterations) + 'iter_best_fitnesses.npy', best_fitness_vals, allow_pickle=True)
        np.save('../results/sat_mut_'+ str(num_seqs)+ 'seq_'+ str(iterations) + 'iter_val_fitnesses.npy', val_fitnesses_by_iter, allow_pickle=True)
        
        np.save('../results/sat_mut_'+ str(num_seqs)+ 'seq_'+ str(iterations) +'iter_seqs.npy', 
                seqs, allow_pickle=True)
    

        t4 = time.time()
        #print(t4-t3)

        print('done')

if __name__ == "__main__":
    main()