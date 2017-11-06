#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Profile-most Probable k-mer Problem
URL: https://stepik.org/lesson/159/step/3?course=Stepic-Interactive-Text-for-Week-3&unit=8217
Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
     Input: A string Text, an integer k, and a 4 × k matrix Profile.
     Output: A Profile-most probable k-mer in Text.
'''

import sys
from functools import reduce

def calculate_max_likelihood_from_dict(inp_str, k, profile):
    'Returns the maximum likelihood k-length sequece within inp_string under the given profile'
    min_str = ''
    max_likelihood = 0.
    for j in range(len(inp_str)-k):
        current_mult = reduce(lambda a,b: a*b, [profile[inp_str[i]][i-j] for i in range(j, j+k)])
        if current_mult > max_likelihood:
            max_likelihood = current_mult
            min_str = inp_str[j:j+k]
    return min_str

if __name__ == '__main__':
    inp_string = sys.stdin.readline()[:-1]
    k = int(sys.stdin.readline())
    rws = [list(map(float, sys.stdin.readline()[:-1].split(" "))) for _ in range(4)]
    keys = ['A','C','G','T']
    profile = dict(zip(keys, rws))
    print(calculate_max_likelihood_from_dict(inp_string, k, profile))

