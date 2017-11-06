#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Implement GreedyMotifSearch
URL: https://stepik.org/lesson/159/step/5?course=Stepic-Interactive-Text-for-Week-3&unit=8217
Code Challenge: Implement GreedyMotifSearch.
    Input: Integers k and t, followed by a collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t).
    If at any step you find more than one Profile-most probable k-mer in a given string, use the
    one occurring first.
'''

import sys
import utils
from functools import reduce


def calculate_seq_likelihood(seq, i_ini, k, profile):
    'Returns the probability asociated to seq[i_ini:i_ini+k] given a profile'
    return reduce(lambda a,b: a*b, [profile[seq[i_ini+i]][i] for i in range(k)])

def get_pos_max_likelihood(inp_str, k, profile):
    'Returns the position of within the input string and the profile of the max likely k-mer'
    min_str_pos = 0
    max_likelihood = 0.
    for j in range(len(inp_str)-k):
        current_mult = calculate_seq_likelihood(inp_str, j, k, profile)
        if current_mult > max_likelihood:
            max_likelihood = current_mult
            min_str_pos = j
    return min_str_pos

def update_counts(current_counts_profile, k, base_pivot, inpstring):
    for i_kmer in range(k):
        current_counts_profile[inpstring[base_pivot+i_kmer], i_kmer] += 1
    return current_counts_profile

def generate_profile(current_counts_profile, pseudocounts = False):
    if pseudocounts:
        aux = current_counts_profile + 1
    else:
        aux = current_counts_profile        
    return np.divide(aux, aux.sum(axis=0))

def GreedyMotifSearch(input_str, k, t, pseudocounts = False):
    'Returns the motifs found in each dna sequence by the greedy method'
    best_seq = []
    best_score = 0

    inputs = [[utils.nucleotid2index[x] for x in seq] for seq in input_str]
    for i_outer in range(len(inputs[0])-k+1):
        current_counts_profile = np.zeros([4, k])
        current_counts_profile = update_counts(current_counts_profile, k, i_outer, inputs[0])
        best_inds = [i_outer]
        for t_seq in range(1, len(inputs)):
            current_profile = generate_profile(current_counts_profile, pseudocounts)
            pos_ind_max = get_pos_max_likelihood(inputs[t_seq], k, current_profile)
            current_counts_profile = update_counts(current_counts_profile, k, pos_ind_max, inputs[t_seq])
            best_inds.append(pos_ind_max)
        # calculate score
        final_profile = np.divide(current_counts_profile, current_counts_profile.sum(axis=0))
        prob = sum(calculate_seq_likelihood(inputs[t_seq], best_inds[t_seq], k, final_profile) for t_seq in range(0, len(inputs)))
        if prob > best_score:
            best_score = prob
            best_seq = best_inds

    return [input_str[t_seq][best_seq[t_seq]:best_seq[t_seq]+k] for t_seq in range(len(inputs))]


if __name__ == '__main__':
    k, t = list(map(int,sys.stdin.readline().split(" ")))
    inp_string = [sys.stdin.readline()[:-1] for i in range(t)]
    pseudocounts = False
    result = GreedyMotifSearch(inp_string, k, t, pseudocounts = pseudocounts)
    print(' '.join(result))

