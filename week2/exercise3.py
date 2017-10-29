from code import HammingDistance
import sys

seq1 = sys.stdin.readline()[:-1]
seq2 = sys.stdin.readline()[:-1]

print(HammingDistance(seq1,seq2))