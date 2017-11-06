#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Implement MedianString
URL: https://stepik.org/lesson/158/step/9?course=Stepic-Interactive-Text-for-Week-3&unit=8216
Code Challenge: Implement MedianString.
     Input: An integer k, followed by a collection of strings Dna.
     Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern. (If there are
     multiple such strings Pattern, then you may return any one.)
'''

import sys
from week3.w3_1A import get_common_keys
from week2.w2_2D import KmersWithAlterations, _kmersWithAlterations_at_index

# input: genome string, k-mer length, max hamming distance allower
# output: set of <kmer, minimum hamming distance at any point of the input sequenceto sequence>
def FindKmersAndHamming(text, k, d):
    'Return set of <kmer, minimum hamming distance at any point of the input sequenceto sequence>'
    min_ham_seqs = {}
    current_seqs = {}
    for i in range(len(text)):
        current_seqs = _kmersWithAlterations_at_index(text, i, current_seqs, k, d)
        length_sequences = len(list(current_seqs.keys())[0])
        if length_sequences == k: # otherwise just iterate (warm-up part of the algorithm)
            for current_seq in current_seqs:
                ham_dist = current_seqs[current_seq]
                if current_seq in min_ham_seqs:
                    min_ham_seqs[current_seq] = min(min_ham_seqs[current_seq], ham_dist)
                else:
                    min_ham_seqs[current_seq] = ham_dist
    return min_ham_seqs


# Implement MedianString.
#    Input: An integer k, followed by a collection of strings Dna.
#    Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern. (If there are
#    multiple such strings Pattern, then you may return any one.)
def MedianString(text_list, k):
    'A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern'
    MAX_D = k
    ham_kmers = [FindKmersAndHamming(text, k, MAX_D) for text in text_list]
    candidate_kmers = get_common_keys(ham_kmers, min_occurence=0)
    min_dist = max([len(text) for text in text_list])
    min_dist_candidate = ''
    for candidate in candidate_kmers:
        now_dist = []
        for seq_kmers in ham_kmers:
            if candidate in seq_kmers:
                now_dist.append(seq_kmers[candidate])
            else:
                now_dist.append(k)
        now_dist = np.sum(now_dist)
        if now_dist < min_dist:
            min_dist = now_dist
            min_dist_candidate = candidate
    return min_dist_candidate   


if __name__ == '__main__':
    k = int(sys.stdin.readline())
    patterns = []
    while(1):
        pattern = sys.stdin.readline()[:-1]
        if len(pattern) == 0:
            break
        else:
            patterns.append(pattern)
    print(MedianString(patterns, k) )

