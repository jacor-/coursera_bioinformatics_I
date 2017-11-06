#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Implement ﻿ApproximatePatternCount
URL: https://stepik.org/lesson/9/step/6?course=Stepic-Interactive-Text-for-Week-2&unit=8224
Code Challenge: Implement ﻿ApproximatePatternCount.
     Input: Strings Pattern and Text as well as an integer d.
     Output: Countd(Text, Pattern).
'''

import sys
import utils
import numpy as np
from week2.w2_2B import ApproximatePatternMatching

def ApproximatePatternCount(seq, pattern, d):
	'Returns the number of times pattern appers in sequence with a Hamming distance <= d'
    patterns = ApproximatePatternMatching(seq, pattern, d)
    return len(patterns)

if __name__ == '__main__':
    pattern = sys.stdin.readline()[:-1]
    seq = sys.stdin.readline()[:-1]
    d = int(sys.stdin.readline()[:-1])
    res = ApproximatePatternCount(seq, pattern, d)
    print(res)
