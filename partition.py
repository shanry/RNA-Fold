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


def count_inside_2(s):
    assert len(s) > 1, "the length of rna should be at least 2!"
    n = len(s)
    counts = np.ones((n, n), dtype=int) # number of structures
    for j in range(1, n):
        for i in range(0, j):
            counts[i, j] = counts[i, j-1]
            if match(s[i], s[j]):
                counts_right = counts[i+1, j-1] if i+1<=j-1 else 1
                for t in range(0, i):
                    counts_left = counts[t, i-1] if t<=i-1 else 1
                    counts[t, j] += counts_left*counts_right
                counts[i, j] += counts_right
    return counts


def count_outside(s, counts=None):
    assert len(s) > 1, "the length of rna should be at least 2!"
    if counts is None:
        counts = count_inside(s)
    assert len(s) == len(counts) and len(s)==len(counts.T), "the length of rna should match counts matrix!"
    n = len(s)
    counts_out = np.zeros((n, n), dtype=int) # number of structures
    counts_out[0, -1] = 1
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


def count_outside_2(s, counts=None):
    assert len(s) > 1, "the length of rna should be at least 2!"
    if counts is None:
        counts = count_inside(s)
    assert len(s) == len(counts) and len(s)==len(counts.T), "the length of rna should match counts matrix!"
    n = len(s)
    counts_out = np.zeros((n, n), dtype=int) # number of structures
    counts_out[0, -1] = 1
    for i in range(1, n):
        counts_out[i, n-1] = counts[0, i-1]
    for j in range(n-1, 0, -1):
        for i in range(0, j):
            counts_out[i, j-1] += counts_out[i, j]
            if match(s[i], s[j]):
                counts_right = counts[i+1, j-1] if i+1<=j-1 else 1
                for t in range(0, i):
                    counts_outleft = counts_out[t, j]
                    counts_out[t, i-1] += counts_outleft*counts_right  # pop left
                counts_outright = counts_out[i, j]
                for t in range(i+2, j):
                    counts_left = counts[i+1, t-1]
                    counts_out[t, j-1] += counts_outright*counts_left # pop right
                if i+1<=j-1:
                    counts_out[i+1, j-1] += counts_outright
    return counts_out


def count_outside_supp(s, counts=None):  # algorithm in the supplementary of linear partition
    assert len(s) > 1, "the length of rna should be at least 2!"
    if counts is None:
        counts = count_inside(s)
    assert len(s) == len(counts) and len(s)==len(counts.T), "the length of rna should match counts matrix!"
    n = len(s)
    counts_out = np.zeros((n, n), dtype=int) # number of structures
    counts_out[0, -1] = 1
    for i in range(1, n):
        counts_out[i, n-1] = counts[0, i-1]
    for j in range(n-1, 0, -1):
        for i in range(0, j):
            counts_out[i, j-1] += counts_out[i, j]
            if match(s[i], s[j]):
                counts_right = counts[i+1, j-1] if i+1<=j-1 else 1
                for t in range(0, i):
                    counts_outleft = counts_out[t, j]
                    counts_out[t, i-1] += counts_outleft*counts_right  # pop left
                    if i+1<=j-1:
                        counts_out[i+1, j-1] += counts_outleft*counts[t, i-1] # pop right
    return counts_out


def pair_match(ss):
    index_pairs = []
    stack = []
    for i, c in enumerate(ss):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            index_pairs.append((j, i))
        else:
            index_pairs.append((i, i))
    return index_pairs

def counts_per_pair(ss_list):
    assert len(ss_list) > 0
    n = len(ss_list[0])
    counts_pp = np.zeros((n, n), dtype=int)
    for i, ss in enumerate(ss_list):
        index_pairs = pair_match(ss)
        for pair in index_pairs:
            counts_pp[pair[0], pair[1]] += 1
    return counts_pp

def verify(p_pair, p_in, p_out):
    n = len(p_pair)
    for i in range(n):
        for j in range(i, n):
            if j==i:
                assert p_pair[i, i] == p_out[i, i] # \beta(i,i) = # of derivations that contain unpaired i
            elif p_pair[i, j] > 0:
                if j==i+1:
                    assert p_pair[i, j] == p_out[i, j] # \beta(i,j) = # of derivations that contain (i,j)
                else:
                    assert p_pair[i, j] == p_in[i+1, j-1]*p_out[i, j] # \alpha(i+1,j-1) * \beta(i,j) = # of derivations that contain (i,j)


def test(algo_in, algo_out):
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
        test_counts_in_out(algo_in, algo_out, rna)
        

def test_counts_in_out(algo_in, algo_out, rna):
    p_in = algo_in(rna)
    print('rna:', rna)
    print('count inside:')
    print(p_in)
    p_out = algo_out(rna, p_in)
    print('count outside:')
    print(p_out)

    best_1, num_struct, best_all = algo_all(rna, None)
    print('best:')
    print(best_1)
    pairs = pair_match(best_1[-1])
    print('mathed pairs:')
    print(pairs)

    counts_pp = counts_per_pair([item[-1] for item in best_all])
    print('count per pair:')
    print(counts_pp)

    verify(counts_pp, p_in, p_out)
    print("verified!")
    
    
def compare_counts(rna):
    p_in = count_inside(rna)
    print('rna:', rna)
    print('count inside:')
    print(p_in)
    
    p_in_2 = count_inside_2(rna)
    print('count inside 2:')
    print(p_in_2)
    assert np.array_equal(p_in, p_in_2)
    
    p_out = count_outside(rna, p_in)
    print('count outside:')
    print(p_out)
    
    p_out_2 = count_outside_2(rna, p_in_2)
    print('count outside 2:')
    print(p_out_2)
    assert np.array_equal(p_out, p_out_2)
    
    print("verified!")
    
    
def compare_counts_supp(rna):
    p_in = count_inside(rna)
    print('rna:', rna)
    print('count inside:')
    print(p_in)
    
    p_in_2 = count_inside_2(rna)
    print('count inside 2:')
    print(p_in_2)
    assert np.array_equal(p_in, p_in_2)
    
    p_out = count_outside(rna, p_in)
    print('count outside:')
    print(p_out)
    
    p_out_supp = count_outside_supp(rna, p_in_2)
    print('count outside supp:')
    print(p_out_supp)
    assert np.array_equal(p_out, p_out_supp)
    
    print("verified!")
    

def count(algo_in, algo_out, rna):
    p_in = algo_in(rna)
    print('rna:', rna)
    print('count inside:')
    print(p_in)
    p_out = algo_out(rna, p_in)
    print('count outside:')
    print(p_out)
    return p_in, p_out


def inside_outside(algo_in, algo_out, rna):
    n = len(rna)
    p_in = algo_in(rna)
    p_out = algo_out(rna, p_in)
    print('count per pair:')
    for i in range(n):
        for j in range(i+1, n):
            if match(rna[i], rna[j]) > 0:
                if j==i+1:
                    count_pair = p_out[i, j] # \beta(i,j) = # of derivations that contain (i,j)
                else:
                    count_pair = p_in[i+1, j-1]*p_out[i, j] # \alpha(i+1,j-1) * \beta(i,j) = # of derivations that contain (i,j)
                print(f"{i+1} {j+1} {count_pair}")
    print('count per unpaired:')
    for i in range(n):
        count_unpair = p_out[i, i]
        print(f"{i+1} {count_unpair}")
    # best_1, num_struct, best_all = algo_all(rna, None)
    # counts_pp = counts_per_pair([item[-1] for item in best_all])
    # print('count per pair:')
    # n = len(rna)
    # for i in range(n):
    #     for j in range(i+1, n):
    #         if counts_pp[i, j] > 0:
    #             print(f"{i+1} {j+1} {counts_pp[i, j]}")
    # print('count per unpaired:')
    # for i in range(n):
    #     print(f"{i+1} {counts_pp[i, i]}")
    # p_in = algo_in(rna)
    # print('rna:', rna)
    # print('count inside:')
    # for i in range(n):
    #     for j in range(i, n):
    #         print(f"{i+1} {j+1} {p_in[i, j]}")
    # p_out = algo_out(rna, p_in)
    # print('count outside:')
    # for i in range(n):
    #     for j in range(i, n):
    #         if p_out[i, j] > 0:
    #             print(f"{i+1} {j+1} {p_out[i, j]}")
    # return p_in, p_out
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", type=str, default="GCACG")
    # parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--algo", type=int, default=1)
    parser.add_argument("--test", action='store_true')
    parser.add_argument("--count", action='store_true')
    parser.add_argument("--inout", action='store_true')
    parser.add_argument("--verify", action='store_true')
    parser.add_argument("--compare", action='store_true')
    args = parser.parse_args()
    print('args:')
    print(args)
    if args.algo == 1:
            partition_in = count_inside
            partition_out = count_outside
    elif args.algo == 2:
        partition_in = count_inside_2
        partition_out = count_outside_2
    else:
        partition_in = count_inside_2
        partition_out = count_outside_supp
    if args.count:
        count(partition_in, partition_out, args.rna)
    if args.inout:
        inside_outside(partition_in, partition_out, args.rna)
    if args.verify:
        test_counts_in_out(partition_in, partition_out, args.rna)
    if args.compare:
        compare_counts(args.rna)
        compare_counts_supp(args.rna)
    if args.test:
        test(partition_in, partition_out)