import argparse

import numpy as np

from k_best import algo_all

np.set_printoptions(linewidth=200)

pairs = {'au', 'gc', 'gu'}

def match(x, y):
    return (x+y).lower() in pairs or (y+x).lower() in pairs


def count_inside(s):
    assert len(s) > 1, "the length of rna should be at least 2!"
    n = len(s)
    counts = np.ones((n, n), dtype=int) # number of structures
    for l in range(1, n):
        for i in range(0, n-l):
            j = i + l
            
            # Case 1 doesn't exist here because l starts from 1
            # Case 2
            counts[i, j] = counts[i, j-1]
            
            # Case 3
            for t in range(i, j):
                if match(s[t],  s[j]):
                        
                    # deal with number of structures
                    count_left = counts[i, t-1] if i<(t-1) else 1
                    count_right = counts[t+1, j-1] if (t+1)<(j-1) else 1
                    counts[i, j] += count_left*count_right
                    
    return counts


def count_outside(s, counts=None):
    assert len(s) > 1, "the length of rna should be at least 2!"
    if counts is None:
        counts = count_inside(s)
    assert len(s) == len(counts) and len(s)==len(counts.T), "the length of rna should match counts matrix!"
    n = len(s)
    counts_out = np.ones((n, n), dtype=int) # number of structures
    for l in range(n-2, -1, -1):
        for i in range(0, n-l):
            j = i + l
            # print(f"l: {l}, i: {i}, j: {j}")
            if j==n-1:
                counts_out[i, j] = counts[0, i-1]
            elif i==0:
                counts_out[i, j] = counts[j+1, n-1]
            else:
                # j+1 unpaired
                counts_out[i, j] = counts_out[i, j+1]
                # print(f'unpaired: {counts[i, j+1]}')
                # j+1 paired rightward
                for t in range(j+2, n):
                    if match(s[j+1],  s[t]):
                        counts_out[i, j] += counts_out[i, t]*counts[j+2, t-1]
                        # print(f'rightward: {counts_out[i, t]*counts[j+2, t-1]}')
                # j+1 paired leftward
                for t in range(0, i):
                    if match(s[t],  s[j+1]):
                        counts_out[i, j] += counts_out[t, j+1]*counts[t+1, i-1]
                        # print(f'leftward: {counts_out[t, j+1]*counts[t+1, i-1]}')
    return counts_out


def pair_match(ss):
    index_pairs = []
    stack = []
    for i, c in enumerate(ss):
        if c == "(":
            stack.append(i)
        if c == ")":
            j = stack.pop()
            index_pairs.append((j, i))
    return index_pairs
    


def test(algo):
    rnas = [
        "ACAGU",
        "AC",
        "GUAC",
        "GCACG",
        "CCGG",
        "CCCGGG",
        "UUCAGGA",
        "AUAACCUA",
        "UUGGACUUG",
        "UUUGGCACUA",
        "GAUGCCGUGUAGUCCAAAGACUUC",
        "AGGCAUCAAACCCUGCAUGGGAGCG"
    ]
    print()
    for i, rna in enumerate(rnas):
        p_in = algo(rna)
        print(i, rna, p_in[0, -1])
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", type=str, default="GCACG")
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--algo", type=int, default=0)
    parser.add_argument("--test", action='store_true')
    args = parser.parse_args()
    print('args:')
    print(args) 
    if args.test:
        test(count_inside)
    else:
        p_in = count_inside(args.rna)
        print('rna:', args.rna)
        print('count inside:')
        print(p_in)
        p_out = count_outside(args.rna, p_in)
        print('count outside:')
        print(p_out)
        
        best_1, num_struct, best_all = algo_all(args.rna, None)
        print('best:')
        print(best_1)
        pairs = pair_match(best_1[-1])
        print('mathed pairs:')
        print(pairs)