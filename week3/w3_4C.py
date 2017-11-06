#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Implement GreedyMotifSearch with pseudocounts
URL: https://stepik.org/lesson/160/step/9?course=Stepic-Interactive-Text-for-Week-3&unit=8218
Code Challenge: Implement GreedyMotifSearch with pseudocounts.
     Input: Integers k and t, followed by a collection of strings Dna.
     Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t) with pseudocounts.
     If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
'''

import sys
import utils
from week3.w3_3B import GreedyMotifSearch

if __name__ == '__main__':
    k, t = list(map(int,sys.stdin.readline().split(" ")))
    inp_string = [sys.stdin.readline()[:-1] for i in range(t)]
    pseudocounts = True
    result = GreedyMotifSearch(inp_string, k, t, pseudocounts = pseudocounts)
    print(' '.join(result))

