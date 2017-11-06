#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: locate CTTGATCAT as a substring in the Vibrio cholerae genome
URL: https://stepik.org/lesson/Some-Hidden-Messages-are-More-Surprising-than-Others-3/step/6?course=Stepic-Interactive-Text-for-Week-1&unit=8232
Code Challenge: Exercise Break (Honors Track): 
	Return a space-separated list of starting positions (in increasing order) where CTTGATCAT appears as a substring in the Vibrio cholerae genome. [Vibrio cholerae Genome (1.1 MB)]

'''

from week1.w1_3B import FindPattern
import sys

if __name__ == '__main__':
    pattern = 'CTTGATCAT'
    sequence = sys.stdin.readline()[:-1]
    print(' '.join(map(str, FindPattern(sequence, pattern))))
