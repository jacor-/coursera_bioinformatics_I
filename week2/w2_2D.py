#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Solve the Frequent Words with Mismatches Problem.
URL: https://stepik.org/lesson/9/step/7?course=Stepic-Interactive-Text-for-Week-2&unit=8224
Code Challenge: Frequent Words with Mismatches Problem: Find the most frequent k-mers with mismatches in a string.
     Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
     Output: All most frequent k-mers with up to d mismatches in Text.
'''

import sys
import utils
import numpy as np


####
## This code is quite convoluted... but I think it is finally quite good. I should clean it anyway :)
## Dynamic programming + recursion
#### 

# input: dictionary structure <sequences length N, hamming_distance from original sequence at index curr_index>
# output: new dictionary <sequences length N-1, hamming_distance from original sequence at index curr_index+1>
# logic: when the first element is drop, the N-1 length resulting sequence can be reconstructed without looking at the whole sequence again
def _update_prev_list(inp, curr_index, prev_combs, length, max_hamming):
    'Returns the list of current possible d-distance sequences by eliminating the previous element of each sequence, merging equal resulting subsequences and updating hamming distances'
    # Propagate result
    if len(list(prev_combs.keys())[0]) != length:
        current_combs = prev_combs
    # First, we clean current_combinations (we drop first character)
    else:
        current_combs = {}
        for now_comb in prev_combs:
            # New combinations drop the first character (slide window)
            new_combination = now_comb[1:]
            if now_comb[0] == inp[max(curr_index-length,0)]: # hamming distance does not change
                new_distance = prev_combs[now_comb]
                current_combs[utils.seq2str(new_combination)] = min(current_combs.setdefault(utils.seq2str(new_combination), max_hamming), new_distance)
            else: # the new sequence has shorter hamming entropy (we drop a difference)
                current_combs[utils.seq2str(new_combination)] = min(current_combs.setdefault(utils.seq2str(new_combination), max_hamming), prev_combs[now_comb]-1)
    return current_combs

# input: dictionary structure <sequences, hamming_distance from original sequence at index curr_index>
# output: new dictionary <sequences, hamming_distance from original sequence at index curr_index+1>
# logic: dinamic programming approach. You dont need to look at the whole sequence, just drop its first element
#        and update the hamming distance until that point
def _kmersWithAlterations_at_index(inp, curr_index, prev_combs, length, max_hamming, alphabet = 'ATCG'):
    'Returns all the kmers with alterations at curr_index. It calculates from previous data, so we do not need to start from scratch at every index position'
    new_combs = {}
    if len(prev_combs) == 0:
        # Initialize
        if max_hamming > 0:
            for alph in alphabet:
                new_combs[alph] = 1
        new_combs[inp[0]] = 0
    else:
        # First, we eliminate the first character of each possible sequence so far
        current_combs = _update_prev_list(inp, curr_index, prev_combs, length, max_hamming)
        # Then, we generate the new batch by adding one character at the end
        for now_comb in current_combs:
            now_ham = current_combs[now_comb]
            if now_ham == max_hamming:
                new_combination = [x for x in now_comb] + [inp[curr_index]] # no permutate in this case
                new_distance = max_hamming
                new_combs[utils.seq2str(new_combination)] = new_combs.setdefault(utils.seq2str(new_combination), max_hamming)
            else:
                for alp in alphabet:
                    new_combination = [x for x in now_comb] + [alp] # no permutate in this case
                    if alp == inp[curr_index]:
                        new_distance = now_ham
                    else:
                        new_distance = now_ham + 1
                    new_combs[utils.seq2str(new_combination)] = min(new_combs.setdefault(utils.seq2str(new_combination), max_hamming), new_distance)
    return new_combs

# input: genome string, k-mer length, max hamming distance allower
# output: count of k-mers with maximum d-hamming distance within the input genome string
def KmersWithAlterations(text, k, d):
    'Returns the number of occurrences of all the k-mers in text with d or more alterations'
    seqcount = {}
    current_seqs = {}
    for i in range(len(text)):
        current_seqs = _kmersWithAlterations_at_index(text, i, current_seqs, k, d)
        if len(list(current_seqs.keys())[0]) == k:
            for current_seq in current_seqs:
                seqcount[''.join(current_seq)] = seqcount.setdefault(''.join(current_seq), 0) + 1
    return seqcount

# input: genome string, k-mer length, max hamming distance allower
# output: most frequent k-mers with maximum d-hamming distance within the input genome string
def FrequentWordsWithMismatches(text, k, d):
    'Returns the most frequent k-mers with at most d mismatches in the input text'
    seqcount = KmersWithAlterations(text, k, d)
    max_counts = np.max(list(seqcount.values()))
    return [x for x in seqcount if seqcount[x] == max_counts]

if __name__ == '__main__':
    inp = sys.stdin.readline()[:-1]
    k, d = map(int, sys.stdin.readline()[:-1].split(' '))
    result = FrequentWordsWithMismatches(inp, k, d)
    print(' '.join(map(str, result)))