#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Find the reverse complement of a DNA string
URL: https://stepik.org/lesson/3/step/2?course=Stepic-Interactive-Text-for-Week-1&unit=8232
Code Challenge: Solve the Pattern Matching Problem.
     Input: Two strings, Pattern and Genome.
     Output: A collection of space-separated integers specifying all starting positions where Pattern appears
     as a substring of Genome.
'''

import numpy as np
from scipy import signal
import utils
import sys

def FindPattern(text, pattern):
    'Returns the starting index where pattern appears as substring of text'
    ref_length = len(pattern)
    text_seq = _convert_to_matrices(text)
    pattern_seq = _convert_to_matrices(pattern)
    matches = signal.correlate(text_seq, pattern_seq, mode = 'valid').flatten() == ref_length
    return [i for i,x in enumerate(matches) if x == 1]

if __name__ == '__main__':
    pattern = sys.stdin.readline()[:-1]
    sequence = sys.stdin.readline()[:-1]
    print(' '.join(map(str, FindPattern(sequence, pattern))))
