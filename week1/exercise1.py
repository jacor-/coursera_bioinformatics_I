from code import PatternCount
import sys

text = sys.stdin.readline()[:-1]
pattern = sys.stdin.readline()[:-1]

print(PatternCount(text, pattern))