import numpy as np
import numpy as np
from scipy import signal
from functools import reduce


seq2str = lambda seq: ''.join(seq)
nucleotid2index = {'G':0,'T':1,'A':2,'C':3}
nucleotid2complement = {'G':'C','C':'G','T':'A','A':'T'}


####
## This code is quite convoluted... but I think it is finally quite good. I should clean it anyway :)
## Dynamic programming + recursion
#### 

# input: dictionary structure <sequences length N, hamming_distance from original sequence at index curr_index>
# output: new dictionary <sequences length N-1, hamming_distance from original sequence at index curr_index+1>
# logic: when the first element is drop, the N-1 length resulting sequence can be reconstructed without looking at the whole sequence again
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


# input: dictionary structure <sequences, hamming_distance from original sequence at index curr_index>
# output: new dictionary <sequences, hamming_distance from original sequence at index curr_index+1>
# logic: dinamic programming approach. You dont need to look at the whole sequence, just drop its first element
#        and update the hamming distance until that point
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

# input: genome string, k-mer length, max hamming distance allower
# output: count of k-mers with maximum d-hamming distance within the input genome string
def FrequencyOfWordsWithAlterations(text, k, d):
    seqcount = {}
    current_seqs = {}
    for i in range(len(text)):
        current_seqs = CalculateKMersAlterations(text, i, current_seqs, k, d)
        if len(list(current_seqs.keys())[0]) == k:
            for current_seq in current_seqs:
                seqcount[''.join(current_seq)] = seqcount.setdefault(''.join(current_seq), 0) + 1
    return seqcount

# input: genome string, k-mer length, max hamming distance allower
# output: most frequent k-mers with maximum d-hamming distance within the input genome string
def FrequentWordsAlterations(text, k, d):
    seqcount = FrequencyOfWordsWithAlterations(text, k, d)
    max_counts = np.max(list(seqcount.values()))
    return [x for x in seqcount if seqcount[x] == max_counts]


# input: a list of dictionaries
# ouput: list of keys appearing in at least min_occurence dictionaries (default min_occurrence = all dictionaries)
def get_common_keys(list_dict, min_occurence = None):
    if min_occurence is None:
        min_occurence = len(list_dict)
    l0 = {}
    for i in range(len(list_dict)):
        li = list_dict[i]
        for key in li.keys():
            if key not in l0:
                l0[key] = 1
            else:
                l0[key] = l0[key] + 1
    return [key for key in l0 if l0[key] >= min_occurence]


# input: genome string, k-mer length, max hamming distance allower
# output: set of <kmer, minimum hamming distance at any point of the input sequenceto sequence>
def FindKmersAndHamming(text, k, d):
    min_ham_seqs = {}
    current_seqs = {}
    for i in range(len(text)):
        current_seqs = CalculateKMersAlterations(text, i, current_seqs, k, d)
        length_sequences = len(list(current_seqs.keys())[0])
        if length_sequences == k: # otherwise just iterate (warm-up part of the algorithm)
            for current_seq in current_seqs:
                ham_dist = current_seqs[current_seq]
                if current_seq in min_ham_seqs:
                    min_ham_seqs[current_seq] = min(min_ham_seqs[current_seq], ham_dist)
                else:
                    min_ham_seqs[current_seq] = ham_dist
    return min_ham_seqs


# Implement MedianString.
#    Input: An integer k, followed by a collection of strings Dna.
#    Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern. (If there are
#    multiple such strings Pattern, then you may return any one.)
def MedianString(text_list, k):
    MAX_D = k
    ham_kmers = [FindKmersAndHamming(text, k, MAX_D) for text in text_list]
    candidate_kmers = get_common_keys(ham_kmers, min_occurence=0)
    min_dist = max([len(text) for text in text_list])
    min_dist_candidate = ''
    for candidate in candidate_kmers:
        now_dist = []
        for seq_kmers in ham_kmers:
            if candidate in seq_kmers:
                now_dist.append(seq_kmers[candidate])
            else:
                now_dist.append(k)
        now_dist = np.sum(now_dist)
        if now_dist < min_dist:
            min_dist = now_dist
            min_dist_candidate = candidate
    return min_dist_candidate   

# maximum likelihood k-length sequece within inp_string under the given profile
def calculate_max_likelihood_from_dict(inp_str, k, profile):
    min_str = ''
    max_likelihood = 0.
    for j in range(len(inp_str)-k):
        current_mult = reduce(lambda a,b: a*b, [profile[inp_str[i]][i-j] for i in range(j, j+k)])
        if current_mult > max_likelihood:
            max_likelihood = current_mult
            min_str = inp_str[j:j+k]
    return min_str


###
# 
###

# probability asociated to seq[i_ini:i_ini+k] given a profile
def calculate_seq_likelihood(seq, i_ini, k, profile):
    return reduce(lambda a,b: a*b, [profile[seq[i_ini+i]][i] for i in range(k)])

# get the position of within the input string and the profile of the max likely k-mer
def get_pos_max_likelihood(inp_str, k, profile):
    min_str_pos = 0
    max_likelihood = 0.
    for j in range(len(inp_str)-k):
        current_mult = calculate_seq_likelihood(inp_str, j, k, profile)
        if current_mult > max_likelihood:
            max_likelihood = current_mult
            min_str_pos = j
    return min_str_pos

def update_counts(current_counts_profile, k, base_pivot, inpstring):
    for i_kmer in range(k):
        current_counts_profile[inpstring[base_pivot+i_kmer], i_kmer] += 1
    return current_counts_profile

def generate_profile(current_counts_profile, pseudocounts = False):
    if pseudocounts:
        aux = current_counts_profile + 1
    else:
        aux = current_counts_profile        
    return np.divide(aux, aux.sum(axis=0))

def GreedyMotifSearch(input_str, k, t, pseudocounts = False):
    best_seq = []
    best_score = 0

    inputs = [[nucleotid2index[x] for x in seq] for seq in input_str]
    for i_outer in range(len(inputs[0])-k+1):
        current_counts_profile = np.zeros([4, k])
        current_counts_profile = update_counts(current_counts_profile, k, i_outer, inputs[0])
        best_inds = [i_outer]
        for t_seq in range(1, len(inputs)):
            current_profile = generate_profile(current_counts_profile, pseudocounts)
            pos_ind_max = get_pos_max_likelihood(inputs[t_seq], k, current_profile)
            current_counts_profile = update_counts(current_counts_profile, k, pos_ind_max, inputs[t_seq])
            best_inds.append(pos_ind_max)
        # calculate score
        final_profile = np.divide(current_counts_profile, current_counts_profile.sum(axis=0))
        prob = sum(calculate_seq_likelihood(inputs[t_seq], best_inds[t_seq], k, final_profile) for t_seq in range(0, len(inputs)))
        if prob > best_score:
            best_score = prob
            best_seq = best_inds

    return [input_str[t_seq][best_seq[t_seq]:best_seq[t_seq]+k] for t_seq in range(len(inputs))]
