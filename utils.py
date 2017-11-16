import numpy as np

nucleotid2index = {'G':0,'T':1,'A':2,'C':3}

def _convert_to_matrices(sequence):
    'Return a matrix <seq, 4> where <i, seq[i]> = 1 and 0 otherwise (indicator of nucleotid at position i)'
    phld = np.zeros((len(sequence), 4))
    for i,s in enumerate(sequence):
        phld[i,utils.nucleotid2index[s]] = 1
    return phld

def HammingDistance(seq1, seq2):
    'Return the Hamming distance between equal-length sequences.'
    if len(seq1) != len(seq2):
        raise ValueError('Undefined for sequences of unequal length.')
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def ReverseComplement(sequence):
    'Returns the reverse complement of a nucleotid sequence'
    return ''.join([nucleotid2complement[sequence[len(sequence)-i-1]] for i in range(len(sequence))])

def seq2str(seq):
    return ''.join(seq)