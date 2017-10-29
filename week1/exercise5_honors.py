from code import FindPattern
import sys

pattern = 'CTTGATCAT'
sequence = sys.stdin.readline()[:-1]

print(' '.join(map(str, FindPattern(sequence, pattern))))
