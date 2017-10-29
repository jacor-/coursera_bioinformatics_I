from code import FrequentWordsAlterations
import sys


inp = sys.stdin.readline()[:-1]
k, d = map(int, sys.stdin.readline()[:-1].split(' '))
result = FrequentWordsAlterations(inp, k, d)

print(' '.join(map(str, result)))