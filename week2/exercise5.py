
from code import ApproximatePatternMatching
import sys

pattern = sys.stdin.readline()[:-1]
seq = sys.stdin.readline()[:-1]
d = int(sys.stdin.readline()[:-1])

res = ApproximatePatternMatching(seq, pattern, d)
print(len(res))