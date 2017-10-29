from code import GreedyMotifSearch
import sys

k, t = list(map(int,sys.stdin.readline().split(" ")))
inp_string = [sys.stdin.readline()[:-1] for i in range(t)]
pseudocounts = False

result = GreedyMotifSearch(inp_string, k, t, pseudocounts = pseudocounts)

print(' '.join(result))
