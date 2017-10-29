from code import FrequentWords
import sys

text = sys.stdin.readline()[:-1]
k = int(sys.stdin.readline()[:-1])

print(' '.join(FrequentWords(text, k)))