from code import FrequentWordsAlterationsWithComplement
import sys

inp = sys.stdin.readline()[:-1]
k, d = map(int, sys.stdin.readline()[:-1].split(' '))

solution = FrequentWordsAlterationsWithComplement(inp, k, d)
print(' '.join(solution))