#!/usr/bin/env python
'''
Content: solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Associated textbook: "Bioinformatics Algorithms: An Active-Learning Approach" by Phillip Compeau & Pavel Pevzner.
Assignment: hosted in Stepik.org
Problem Title: Find the reverse complement of a DNA string
URL: https://stepik.org/lesson/3/step/5?course=Stepic-Interactive-Text-for-Week-1&unit=8232
Code Challenge: Reverse Complement Problem: Find the reverse complement of a DNA string.
     Input: A DNA string Pattern.
     Output: Patternrc , the reverse complement of Pattern.
'''

import numpy as np
from scipy import signal
import utils
import sys


if __name__ == '__main__':
    text = sys.stdin.readline()[:-1]
    print(utils.ReverseComplement(text))