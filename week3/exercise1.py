from code import FrequencyOfWordsWithAlterations, get_common_keys
import sys

#Code Challenge: Implement MotifEnumeration (reproduced below).
#Input: Integers k and d, followed by a collection of strings Dna.
#Output: All (k, d)-motifs in Dna.
def MotifEnumeration(Dna, k, d):
	seqs = [FrequencyOfWordsWithAlterations(sequence, k, d) for sequence in Dna]
	return get_common_keys(seqs)


k,d = map(int, sys.stdin.readline()[:-1].split(" "))
patterns = []
while(1):
	pattern = sys.stdin.readline()[:-1]
	if len(pattern) == 0:
		break
	else:
		patterns.append(pattern)


print(' '.join(MotifEnumeration(patterns, k, d) ))