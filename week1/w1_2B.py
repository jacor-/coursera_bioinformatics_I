#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Solve the Frequent Words Problem
URL: https://stepik.org/lesson/2/step/10?course=Stepic-Interactive-Text-for-Week-1&unit=8231
Code Challenge: Solve the Frequent Words Problem.
     Input: A string Text and an integer k.
     Output: All most frequent k-mers in Text.
'''

import numpy as np
from scipy import signal
import utils
import sys


def FrequentWords(text, k):
    'Returns most frequent k-mers in text'
    seqs = {}
    current_seq = [text[i] for i in range(k)]
    seqs = {''.join(current_seq): 1}
    for i in range(k, len(text)):
        current_seq = current_seq[1:] + [text[i]]
        seqs[''.join(current_seq)] = seqs.setdefault(''.join(current_seq), 0) + 1
    max_counts = np.max(list(seqs.values()))
    return [x for x in seqs if seqs[x] == max_counts]

if __name__ == '__main__':
    text = sys.stdin.readline()[:-1]
    k = int(sys.stdin.readline()[:-1])
    print(' '.join(FrequentWords(text, k)))