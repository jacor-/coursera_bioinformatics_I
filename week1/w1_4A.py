#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Clump Finding Problem
URL: https://stepik.org/lesson/4/step/5?course=Stepic-Interactive-Text-for-Week-1&unit=8233
Code Challenge: Clump Finding Problem: Find patterns forming clumps in a string.
     Input: A string Genome, and integers k, L, and t.
     Output: All distinct k-mers forming (L, t)-clumps in Genome.
'''

import numpy as np
from scipy import signal
import utils
import sys

def _FindPos_AllKMers(sequence,k):
    'Returns a dictionary {kmer: [appearing_positions]} containing all the kmers in sequence'
    seqs = {}
    current_seq = [sequence[i] for i in range(k)]
    seqs = {utils.seq2str(current_seq): [0]}
    for i in range(1, len(sequence)-k+1):
        current_seq = current_seq[1:] + [sequence[i+k-1]]
        if utils.seq2str(current_seq) in seqs.keys():
            seqs[utils.seq2str(current_seq)].append(i)
        else:
            seqs[utils.seq2str(current_seq)] = [i]
    return seqs
    
def _is_sequence_clump(pos_list, L, t, k):
    'Returns if the subsequences starting at pos_list appears t or mode times in a substring of length L'
    for i in range(t-1, len(pos_list)):
        if pos_list[i] - pos_list[i-t+1] < L-k+1:
            return True
    return False

def FindClumps(seq, k, L, t):
    'Returns a list of kmers appearing more than t times in a L-length substring of seq'
    all_kmers = _FindPos_AllKMers(seq, k)
    return [kmer for kmer in all_kmers.keys() if _is_sequence_clump(all_kmers[kmer], L, t, k)]

if __name__ == '__main__':
    sequence = sys.stdin.readline()[:-1]
    k, L, t = map(int, sys.stdin.readline()[:-1].split(' '))
    solution = FindClumps(sequence, k, L, t)
    print(' '.join(solution))