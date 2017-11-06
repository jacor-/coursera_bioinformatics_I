#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Solve the Frequent Words with Mismatches Problem and reverse complement.
URL: https://stepik.org/lesson/9/step/8?course=Stepic-Interactive-Text-for-Week-2&unit=8224
Code Challenge: Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent k-mers (with mismatches and reverse complements) in a string.
      Input: A DNA string Text as well as integers k and d.
      Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Patternrc) over all possible k-mers.
'''

import sys
import utils
import numpy as np
from week2.w2_2D import FrequentWordsWithMismatches

def FrequentWordsMismatchesWithComplement(text, k, d):
    'Returns the most frequent k-mers and reverse complement sequences with at most d mismatches in the input text'
    seqcount = FrequentWordsWithMismatches(text, k, d)
    rev_sum = {}
    for element in seqcount:
        reverse_element = utils.ReverseComplement(element)
        new_value = seqcount[element] + (seqcount[reverse_element] if reverse_element in seqcount else 0)
        rev_sum[element] = new_value
        rev_sum[reverse_element] = new_value
    max_counts = np.max(list(rev_sum.values()))
    return [x for x in rev_sum if rev_sum[x] == max_counts]

if __name__ == '__main__':
	inp = sys.stdin.readline()[:-1]
	k, d = map(int, sys.stdin.readline()[:-1].split(' '))
	solution = FrequentWordsMismatchesWithComplement(inp, k, d)
	print(' '.join(solution))