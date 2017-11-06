#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Find 9-mers forming (500,3)-clumps in E-coli genome
URL: https://stepik.org/lesson/4/step/6?unit=8233
Code Challenge: Exercise Break: How many different 9-mers form (500,3)-clumps in the E. coli genome? (In other words, do not count a 9-mer more than once.)
'''

from week1.w1_4A import FindClumps
import sys
import settings

if __name__ == '__main__':
    ecoli = open(settings.dataset_paths + '/E_coli.txt').readline()
    k, L, t = 9, 500, 3
    solution = FindClumps(ecoli, k, L, t)
    print(len(solution))
