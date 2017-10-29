from code import SkewGenome
import sys

seq = sys.stdin.readline()[:-1]
print(' '.join(map(str,SkewGenome(seq))))