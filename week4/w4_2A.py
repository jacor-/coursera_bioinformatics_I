#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Implement GibbsSamplerMotifs
URL: https://stepik.org/lesson/163/step/4?course=Stepic-Interactive-Text-for-Week-4&unit=8210
Code Challenge: Implement GibbsSamplerMotifs.
     Input: Integers k, t, and N, followed by a collection of strings Dna.
     Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with
     20 random starts. Remember to use pseudocounts!
'''

from random import randint, choice
from utils import HammingDistance
from week4.w4_1A import randomized_motif_search, score, profile_with_pseudocounts
import sys
import numpy

def profile_randomized_kmer(dna, k, prof):
    '''Return the profile most probable k-mer in a given dna sequence.'''
    # A dictionary relating nucleotides to their position within the profile.
    nuc_loc = {nucleotide:index for index,nucleotide in enumerate('ACGT')}
    # Initialize the maximum probabily.
    probs = []
    # Compute the probability of the each k-mer, store it if it's currently a maximum.
    for i in range(len(dna)-k):
        current_prob = 1.
        for j, nucleotide in enumerate(dna[i:i+k]):
            current_prob *= prof[j][nuc_loc[nucleotide]]
        probs.append(current_prob)

    # weighted probs
    i = numpy.random.choice(len(probs), p=numpy.array(probs) / numpy.sum(probs))
    return dna[i:i+k]


def gibbs_sampling_motif_search(dna_list,k,t,N,init_motifs = None):
    if init_motifs:
        motifs = init_motifs
    else:
        # Randomly generate k-mers from each sequence in the dna list. Do N steps in the gibbs sampler
        rand_ints = [randint(0,len(dna_list[0])-k) for a in range(t)]
        motifs = [dna_list[i][r:r+k] for i,r in enumerate(rand_ints)]

    # Initialize the best score as a score higher than the highest possible score.
    best_score = [score(motifs), list(motifs)]

    # Do N steps
    for j in range(N):
        # select the ith motif to be changed
        i = randint(0,t-1)
        # calculate the new profile ignoring the ith motif
        current_profile = profile_with_pseudocounts([x for amotif,x in enumerate(motifs) if amotif != i])
        # find the best motif given that profile and the ith dna string
        motifs[i] = profile_randomized_kmer(dna_list[i], k, current_profile)
        # keep the newly generated motif set if the score improves
        current_score = score(motifs)
        if current_score < best_score[0]:
            best_score = [current_score, list(motifs)]

    return best_score

if __name__ == '__main__':
    
    k,t,N = map(int, sys.stdin.readline().split())
    dna_list = [line.strip() for line in sys.stdin.readlines()]

    # Initialize the best scoring motifs as a score higher than the highest possible score.
    best_motifs = [k*t, None]

    # Then we run the Gibbs sampling chain
    for repeat in range(20):
        current_motifs = gibbs_sampling_motif_search(dna_list,k,t,N)
        if current_motifs[0] < best_motifs[0]:
            best_motifs = current_motifs

    # Print and save the answer.
    print('\n'.join(best_motifs[1]))
    