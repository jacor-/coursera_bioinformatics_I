#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Minimum Skew Problem
URL: https://stepik.org/lesson/7/step/6?course=Stepic-Interactive-Text-for-Week-2&unit=8223
Code Challenge:  Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.
     Input: A DNA string Genome.
     Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).s
'''

import sys
import settings
import numpy as np
from week2.w2_1A import SkewGenome

def FindMinimumSkew(dna):
    'Return the positions of minimum skew given the input dna'
    els = SkewGenome(dna)
    mn = np.min(els)
    return [i for i in range(len(els)) if els[i] == mn]

if __name__ == '__main__':
    dna = sys.stdin.readline()[:-1]
    print(' '.join(map(str,FindMinimumSkew(dna))))