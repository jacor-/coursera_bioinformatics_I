#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Implement MotifEnumeration
URL: https://stepik.org/lesson/156/step/8?course=Stepic-Interactive-Text-for-Week-3&unit=8214
Code Challenge: Implement MotifEnumeration (reproduced below).
     Input: Integers k and d, followed by a collection of strings Dna.
     Output: All (k, d)-motifs in Dna.
'''

import sys
from week2.w2_2E import KmersWithAlterations

# input: a list of dictionaries
# ouput: list of keys appearing in at least min_occurence dictionaries (default min_occurrence = all dictionaries)
def get_common_keys(list_dict, min_occurence = None):
    'Returns keys appearing in at least min_occurence input dictionaries'
    if min_occurence is None:
        min_occurence = len(list_dict)
    l0 = {}
    for i in range(len(list_dict)):
        li = list_dict[i]
        for key in li.keys():
            if key not in l0:
                l0[key] = 1
            else:
                l0[key] = l0[key] + 1
    return [key for key in l0 if l0[key] >= min_occurence]


def MotifEnumeration(Dna, k, d):
    seqs = [KmersWithAlterations(sequence, k, d) for sequence in Dna]
    return get_common_keys(seqs)

if __name__ == '__main__':
    k,d = map(int, sys.stdin.readline()[:-1].split(" "))
    patterns = []
    while(1):
        pattern = sys.stdin.readline()[:-1]
        if len(pattern) == 0:
            break
        else:
            patterns.append(pattern)

    print(' '.join(MotifEnumeration(patterns, k, d) ))


