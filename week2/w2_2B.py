#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Solve the Approximate Pattern Matching Problem
URL: https://stepik.org/lesson/9/step/4?course=Stepic-Interactive-Text-for-Week-2&unit=8224
Code Challenge:      Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
     Input: Strings Pattern and Text along with an integer d.
     Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.


'''

import sys
import utils
import numpy as np
from scipy import signal


def ApproximatePatternMatching(seq, pattern, d):
	'Returns positions in seq where pattern appears with a hamming distance <= d'
    ref_length = len(pattern)
    text_seq = utils._convert_to_matrices(seq)
    pattern_seq = utils._convert_to_matrices(pattern)
    matches = signal.correlate(text_seq, pattern_seq, mode = 'valid').flatten()
    return [i for i,x in enumerate(matches) if x >= ref_length-d]

if __name__ == '__main__':
	pattern = sys.stdin.readline()[:-1]
	seq = sys.stdin.readline()[:-1]
	d = int(sys.stdin.readline()[:-1])
	res = ApproximatePatternMatching(seq, pattern, d)
	print(' '.join(map(str, res)))

