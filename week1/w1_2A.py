#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Implement PatternCount
URL: https://stepik.org/lesson/2/step/7?course=Stepic-Interactive-Text-for-Week-1&unit=8231
Code Challenge: Implement PatternCount (reproduced below).
     Input: Strings Text and Pattern.
     Output: Count(Text, Pattern).
'''

import numpy as np
from scipy import signal
import utils
import sys


def PatternCount(text, pattern):
    ref_length = len(pattern)
    text_seq = utils._convert_to_matrices(text)
    pattern_seq = utils._convert_to_matrices(pattern)
    return np.sum(signal.correlate(text_seq, pattern_seq, mode = 'valid').flatten() == ref_length) 

if __name__ == '__main__':
    text = sys.stdin.readline()[:-1]
    pattern = sys.stdin.readline()[:-1]
    print(PatternCount(text, pattern))