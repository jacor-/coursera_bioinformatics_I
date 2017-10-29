import numpy as np
import numpy as np
from scipy import signal

seq2str = lambda seq: ''.join(seq)
nucleotid2index = {'G':0,'T':1,'A':2,'C':3}
nucleotid2complement = {'G':'C','C':'G','T':'A','A':'T'}

def _convert_to_matrices(sequence):
    phld = np.zeros((len(sequence), 4))
    for i,s in enumerate(sequence):
        phld[i,nucleotid2index[s]] = 1
    return phld

def SkewGenome(seq):
    count = [0]
    for i in range(len(seq)):
        if seq[i] == 'G':
            count.append(count[-1]+1)
        elif seq[i] == 'C':
            count.append(count[-1]-1)
        else:
            count.append(count[-1])
    return count

def FindMinimumSkew(seq):
	els = SkewGenome(seq)
	mn = np.min(els)
	return [i for i in range(len(els)) if els[i] == mn]

def HammingDistance(seq1, seq2):
	dst = 0
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			dst += 1
	return dst

def ApproximatePatternMatching(seq, pattern, d):
    ref_length = len(pattern)
    text_seq = _convert_to_matrices(seq)
    pattern_seq = _convert_to_matrices(pattern)
    matches = signal.correlate(text_seq, pattern_seq, mode = 'valid').flatten()
    return [i for i,x in enumerate(matches) if x >= ref_length-d]


####
## This code is quite convoluted... but I think it is finally quite good. I should clean it anyway :)
## Dynamic programming + recursion
#### 


def _update_prev_list(inp, curr_index, prev_combs, length, max_hamming):
    # Propagate result
    if len(list(prev_combs.keys())[0]) != length:
        current_combs = prev_combs
    # First, we clean current_combinations (we drop first character)
    else:
        current_combs = {}
        for now_comb in prev_combs:
            # New combinations drop the first character (slide window)
            new_combination = now_comb[1:]
            if now_comb[0] == inp[max(curr_index-length,0)]: # hamming distance does not change
                new_distance = prev_combs[now_comb]
                current_combs[seq2str(new_combination)] = min(current_combs.setdefault(seq2str(new_combination), max_hamming), new_distance)
            else: # the new sequence has shorter hamming entropy (we drop a difference)
                current_combs[seq2str(new_combination)] = min(current_combs.setdefault(seq2str(new_combination), max_hamming), prev_combs[now_comb]-1)
    return current_combs

def CalculateKMersAlterations(inp, curr_index, prev_combs, length, max_hamming, alphabet = 'ATCG'):
    new_combs = {}
    if len(prev_combs) == 0:
        # Initialize
        if max_hamming > 0:
            for alph in alphabet:
                new_combs[alph] = 1
        new_combs[inp[0]] = 0
    else:
        current_combs = _update_prev_list(inp, curr_index, prev_combs, length, max_hamming)
        # Then, we generate the new batch by adding one character at the end
        for now_comb in current_combs:
            now_ham = current_combs[now_comb]
            if now_ham == max_hamming:
                new_combination = [x for x in now_comb] + [inp[curr_index]] # no permutate in this case
                new_distance = max_hamming
                new_combs[seq2str(new_combination)] = new_combs.setdefault(seq2str(new_combination), max_hamming)
            else:
                # We can generate caos here...
                for alp in alphabet:
                    new_combination = [x for x in now_comb] + [alp] # no permutate in this case
                    if alp == inp[curr_index]:
                        new_distance = now_ham
                    else:
                        new_distance = now_ham + 1
                    new_combs[seq2str(new_combination)] = min(new_combs.setdefault(seq2str(new_combination), max_hamming), new_distance)
    return new_combs

def FrequencyOfWordsWithAlterations(text, k, d):
    seqcount = {}
    current_seqs = {}
    for i in range(len(text)):
        current_seqs = CalculateKMersAlterations(text, i, current_seqs, k, d)
        if len(list(current_seqs.keys())[0]) == k:
            for current_seq in current_seqs:
                seqcount[''.join(current_seq)] = seqcount.setdefault(''.join(current_seq), 0) + 1
    return seqcount

def FrequentWordsAlterations(text, k, d):
    seqcount = FrequencyOfWordsWithAlterations(text, k, d)
    max_counts = np.max(list(seqcount.values()))
    return [x for x in seqcount if seqcount[x] == max_counts]

def ReverseComplement(sequence):
    return ''.join([nucleotid2complement[sequence[len(sequence)-i-1]] for i in range(len(sequence))])

def FrequentWordsAlterationsWithComplement(text, k, d):
    seqcount = FrequencyOfWordsWithAlterations(text, k, d)

    rev_sum = {}
    for element in seqcount:
        reverse_element = ReverseComplement(element)
        new_value = seqcount[element] + (seqcount[reverse_element] if reverse_element in seqcount else 0)

        rev_sum[element] = new_value
        rev_sum[reverse_element] = new_value

    max_counts = np.max(list(rev_sum.values()))
    return [x for x in rev_sum if rev_sum[x] == max_counts]
