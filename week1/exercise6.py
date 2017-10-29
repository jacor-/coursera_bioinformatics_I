
from code import FindClumps
import sys

sequence = sys.stdin.readline()[:-1]
k, L, t = map(int, sys.stdin.readline()[:-1].split(' '))

solution = FindClumps(sequence, k, L, t)

print(len(solution))
#print(' '.join(solution))