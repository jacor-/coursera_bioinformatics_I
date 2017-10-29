from code import FindMinimumSkew
import sys

seq = sys.stdin.readline()[:-1]
print(' '.join(map(str,FindMinimumSkew(seq))))