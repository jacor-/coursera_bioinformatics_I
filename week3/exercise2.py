from code import MedianString
import sys

#Code Challenge: Implement MotifEnumeration (reproduced below).
#Input: Integers k and d, followed by a collection of strings Dna.
#Output: All (k, d)-motifs in Dna.


k = int(sys.stdin.readline())
patterns = []
while(1):
	pattern = sys.stdin.readline()[:-1]
	if len(pattern) == 0:
		break
	else:
		patterns.append(pattern)


print(MedianString(patterns, k) )