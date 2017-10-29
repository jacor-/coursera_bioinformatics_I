import numpy as np
from scipy import signal


nucleotid2index = {'G':0,'T':1,'A':2,'C':3}
nucleotid2complement = {'G':'C','C':'G','T':'A','A':'T'}

def _convert_to_matrices(sequence):
    phld = np.zeros((len(sequence), 4))
    for i,s in enumerate(sequence):
        phld[i,nucleotid2index[s]] = 1
    return phld

def PatternCount(text, pattern):
    ref_length = len(pattern)
    text_seq = _convert_to_matrices(text)
    pattern_seq = _convert_to_matrices(pattern)

    return np.sum(signal.correlate(text_seq, pattern_seq, mode = 'valid').flatten() == ref_length) 

def FindPattern(text, pattern):
    ref_length = len(pattern)
    text_seq = _convert_to_matrices(text)
    pattern_seq = _convert_to_matrices(pattern)
    matches = signal.correlate(text_seq, pattern_seq, mode = 'valid').flatten() == ref_length
    return [i for i,x in enumerate(matches) if x == 1]

def FrequentWords(text, k):
    seqs = {}
    current_seq = [text[i] for i in range(k)]
    seqs = {''.join(current_seq): 1}
    for i in range(k, len(text)):
        current_seq = current_seq[1:] + [text[i]]
        seqs[''.join(current_seq)] = seqs.setdefault(''.join(current_seq), 0) + 1

    max_counts = np.max(list(seqs.values()))
    return [x for x in seqs if seqs[x] == max_counts]

def ReverseComplement(sequence):
    return ''.join([nucleotid2complement[sequence[len(sequence)-i-1]] for i in range(len(sequence))])


def _FindPos_AllKMers(sequence,k):
    #Let's find all the sequences and the index where they start
    seq2str = lambda seq: ''.join(seq)
    seqs = {}

    current_seq = [sequence[i] for i in range(k)]
    seqs = {seq2str(current_seq): [0]}
    for i in range(1, len(sequence)-k+1):
        current_seq = current_seq[1:] + [sequence[i+k-1]]
        if seq2str(current_seq) in seqs.keys():
            seqs[seq2str(current_seq)].append(i)
        else:
            seqs[seq2str(current_seq)] = [i]
    return seqs
    
def _is_sequence_clump_v2(pos_list, L, t, k):
    # To check if yes or no, we just need to slide a window of length t!
    for i in range(t-1, len(pos_list)):
        if pos_list[i] - pos_list[i-t+1] < L-k+1:
            return True
    return False

def FindClumps(seq, k, L, t):
    all_kmers = _FindPos_AllKMers(seq, k)
    return [kmer for kmer in all_kmers.keys() if _is_sequence_clump_v2(all_kmers[kmer], L, t, k)]
  

    