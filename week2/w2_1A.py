#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Find the skew sequence of a given input
URL: https://stepik.org/lesson/7/step/4?course=Stepic-Interactive-Text-for-Week-2&unit=8223
Code Challenge: Exercise Break: Give all values of Skewi (GAGCCACCGCGATA) for i ranging from 0 to 14.
'''

import sys
import settings

def SkewGenome(dna):
	'Returns the skew sequence based on input dna string'
    count = [0]
    for i in range(len(dna)):
        if dna[i] == 'G':
            count.append(count[-1]+1)
        elif dna[i] == 'C':
            count.append(count[-1]-1)
        else:
            count.append(count[-1])
    return count

if __name__ == '__main__':
	dna = sys.stdin.readline()[:-1]
	print(' '.join(map(str,SkewGenome(dna))))