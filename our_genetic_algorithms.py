from __future__ import print_function, division
import string
import numpy as np
import seq_functions

import tensorflow as tf
from tensorflow.keras.models import load_model
from sklearn import metrics
from Bio import SeqIO
import math
import os
import random

# make it like this
#https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html#sklearn.linear_model.LogisticRegression

def main():
    
    random.seed(100)
    
    # this must be the same order as the one hot encoding!
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    
    #base_rates = [0.2, 0.3, 0.3, 0.2]
    base_rates = [0.1, 0.9, 0, 0]
    assert np.sum(np.asarray(base_rates)) == 1.0
    
    # old algorithm
    #target_string = "Genetic Algorithm"
    #population_size = 100
    
    # new algorithm
    mutation_rate = 1
    crossover_rate = 0
    crossover_length = 5
    
    ## test sequences
    test_seqs = ['AGGATTAGGA', 
                 'TTTTTTTTTT']
    test_seq_list = [[j for j in s] for s in test_seqs]
    print(test_seq_list)
    seqs = np.array([seq_functions.onehot_seq(seq) for seq in test_seq_list])
    
    # load ML model
    model_name = '/home/mleone2/ml_class_project/pretrained_models/Mo2015_EXCpos_Ctx_fold1/Mo2015_EXCpos_Ctx_fold1_OCP_NB1000_NE23_BR0.01_MR0.1_BM0.85_MM0.99_DO0.1.h5'
    model = load_model(model_name, compile=False)
    
    #ga = GeneticAlgorithm_old(target_string,
    #                                    population_size,
    #                                    mutation_rate)
    
    ga = GeneticAlgorithm(mutation_rate, crossover_rate, 
                          crossover_length, base_rates)

    ga.run(iterations=10, sequences = seqs, model = model_name)

    
class GeneticAlgorithm(): 
    def __init__(self, mutation_rate, crossover_rate, crossover_length, base_rates):
        
        
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.crossover_length = crossover_length
        
        self.base_rates = base_rates # list assumed to correspond to order of one hot encoding
        
    def _mutate(self, sequences):
        
        """ Randomly change the individual's characters with probability
        self.mutation_rate """
        
        ### get N X length of sequence binary random vars
        ### assume sequences is shape (N,length of sequence, 4)
        ### 1 at coordinate (x,y) means mutate sequence x at base y
        mutate_yesno = np.random.choice([0, 1], size=sequences.shape[0:2], p=[1-self.mutation_rate, self.mutation_rate])
        
        ### which_base is the set of bases to mutate to
        ### in future if slow: only calculate for bases that will mutate at all
        #which_base = np.random.choice(4, mutate_yesno[mutate_yesno == 1].shape, p=self.base_rates)
        which_base = np.random.choice(4, mutate_yesno.shape, p=self.base_rates)
        
        # generate new sequence after mutating
        new_sequences = np.zeros(sequences.shape)
        # keep non-mutated one hot encoding
        new_sequences[mutate_yesno == 0,:] = sequences[mutate_yesno == 0,:]
        
        # indices where to mutate
        nonzero_rows_columns = np.nonzero(mutate_yesno)
        
        # get new one hot encoding where mutated
        row_iter = 0
        col_iter = 0
        if len(nonzero_rows_columns[0]) > 0:
            for ii in nonzero_rows_columns[0]:
                for jj in nonzero_rows_columns[1]:
                    
                    # in future if slow, use different which_base
                    #new_sequences[ii,jj,
                    #              which_base[row_iter,col_iter]] = 1
                    new_base_id = which_base[ii,jj]
                    new_sequences[ii,jj,new_base_id] = 1
                    
                    col_iter = col_iter + 1
                        
                row_iter = row_iter + 1
                
        return new_sequences
                
    
    def _unonehot(self, sequences):
        ## essentially for testing
        
        letters = ['A','C','G','T']
        string_sequences = []
        for seq_id in range(sequences.shape[0]):
            #seq_onehot = sequences[seq_id,:,:]
            
            new_string = []
            for letter_id in range(sequences.shape[1]):
                which_letter = np.nonzero(sequences[seq_id,
                                                    letter_id,])[0][0]

                new_string.append(letters[which_letter])
                
            string_sequences.append(new_string)
            
        return ["".join(s) for s in string_sequences]
    
    def _predict(self, sequences, model):
        y_pred_score = model.predict(seqs)
            
        
    def run(self, iterations,sequences, model):
        
        print(self._unonehot(sequences))
        for ii in range(iterations):
            
            if multiplier > 1:
                self._mutate(sequences)
                
                
            candidates = self._mutate(sequences)
            #sequences = self._crossover(sequences)
            
            
            # get fitness (model predictions)
            fitness =  model.predict(candidates)
            
            
            print(self._unonehot(sequences))
            
        return sequences
        
        
    
class GeneticAlgorithm_old():
    """An implementation of a Genetic Algorithm which will try to produce the user
    specified target string.
    Parameters:
    -----------
    target_string: string
        The string which the GA should try to produce.
    population_size: int
        The number of individuals (possible solutions) in the population.
    mutation_rate: float
        The rate (or probability) of which the alleles (chars in this case) should be
        randomly changed.
    """
    def __init__(self, target_string, population_size, mutation_rate):
        self.target = target_string
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.letters = [" "] + list(string.ascii_letters)

    def _initialize(self):
        """ Initialize population with random strings """
        self.population = []
        for _ in range(self.population_size):
            # Select random letters as new individual
            individual = "".join(np.random.choice(self.letters, size=len(self.target)))
            self.population.append(individual)

    def _calculate_fitness(self):
        """ Calculates the fitness of each individual in the population """
        population_fitness = []
        for individual in self.population:
            # Calculate loss as the alphabetical distance between
            # the characters in the individual and the target string
            loss = 0
            for i in range(len(individual)):
                letter_i1 = self.letters.index(individual[i])
                letter_i2 = self.letters.index(self.target[i])
                loss += abs(letter_i1 - letter_i2)
            fitness = 1 / (loss + 1e-6)
            population_fitness.append(fitness)
        return population_fitness

    def _mutate(self, individual):
        """ Randomly change the individual's characters with probability
        self.mutation_rate """
        individual = list(individual)
        for j in range(len(individual)):
            # Make change with probability mutation_rate
            if np.random.random() < self.mutation_rate:
                individual[j] = np.random.choice(self.letters)
        # Return mutated individual as string
        return "".join(individual)

    def _crossover(self, parent1, parent2):
        """ Create children from parents by crossover """
        # Select random crossover point
        cross_i = np.random.randint(0, len(parent1))
        child1 = parent1[:cross_i] + parent2[cross_i:]
        child2 = parent2[:cross_i] + parent1[cross_i:]
        return child1, child2

    def run(self, iterations):
        # Initialize new population
        self._initialize()

        for epoch in range(iterations):
            population_fitness = self._calculate_fitness()

            fittest_individual = self.population[np.argmax(population_fitness)]
            highest_fitness = max(population_fitness)

            # If we have found individual which matches the target => Done
            if fittest_individual == self.target:
                break

            # Set the probability that the individual should be selected as a parent
            # proportionate to the individual's fitness.
            parent_probabilities = [fitness / sum(population_fitness) for fitness in population_fitness]

            # Determine the next generation
            new_population = []
            for i in np.arange(0, self.population_size, 2):
                # Select two parents randomly according to probabilities
                parent1, parent2 = np.random.choice(self.population, size=2, p=parent_probabilities, replace=False)
                # Perform crossover to produce offspring
                child1, child2 = self._crossover(parent1, parent2)
                # Save mutated offspring for next generation
                new_population += [self._mutate(child1), self._mutate(child2)]

            print ("[%d Closest Candidate: '%s', Fitness: %.2f]" % (epoch, fittest_individual, highest_fitness))
            self.population = new_population

        print ("[%d Answer: '%s']" % (epoch, fittest_individual))
        
        
if __name__ == "__main__":
    main()