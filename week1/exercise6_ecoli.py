
from code import FindClumps
import sys

ecoli = open('/home/jose/Downloads/E_coli.txt').readline()
k, L, t = 9, 500, 3

solution = FindClumps(ecoli, k, L, t)

print(len(solution))
