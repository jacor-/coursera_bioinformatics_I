#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Hamming Distance Problem
URL: https://stepik.org/lesson/9/step/3?course=Stepic-Interactive-Text-for-Week-2&unit=8224
Code Challenge:      Hamming Distance Problem: Compute the Hamming distance between two strings.
     Input: Two strings of equal length.
     Output: The Hamming distance between these strings.
'''

import sys
import utils
import numpy as np

if __name__ == '__main__':
    seq1 = sys.stdin.readline()[:-1]
    seq2 = sys.stdin.readline()[:-1]
    print(utils.HammingDistance(seq1,seq2))

