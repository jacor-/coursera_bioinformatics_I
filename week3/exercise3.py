from code import calculate_max_likelihood_from_dict
import sys

inp_string = sys.stdin.readline()[:-1]
k = int(sys.stdin.readline())
rws = [list(map(float, sys.stdin.readline()[:-1].split(" "))) for _ in range(4)]
keys = ['A','C','G','T']
profile = dict(zip(keys, rws))

print(calculate_max_likelihood_from_dict(inp_string, k, profile))

