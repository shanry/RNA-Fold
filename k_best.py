import numpy as np

import heapq

pairs = {'au', 'gc', 'gu'}

def match(x, y):
    return (x+y).lower() in pairs or (y+x).lower() in pairs


class MinHeap:
 
    def __init__(self, k):
        self.pool = []
        self.k = k
 
    def add(self, n):
 
        # if the min-heap's size is less than `k`, push directly
        if len(self.pool) < self.k:
            heapq.heappush(self.pool, n)
 
        # otherwise, compare to decide push or not
        elif self.pool[0][0] < n[0]:
            heapq.heappushpop(self.pool, n)
            

def algo_0(s, k=10):
    assert len(s) > 1, "the length of rna should be at least 2!"
    n = len(s)
    counts = np.ones((n, n)) # number of structures
    opts_k = np.zeros((n, n) ,dtype=object)
    for row in range(n):
        for col in range(n):
            opts_k[row, col] = []
            if row == col:
                opts_k[row, col].append((0, '.'))  # tuple: (num_pair, prediction)
    for l in range(1, n):
        for i in range(0, n-l):
            j = i + l
            
            # Case 1 doesn't exist here because l starts from 1
            # Case 2
            counts[i, j] = counts[i, j-1]
            minhp = MinHeap(k)
            for element in opts_k[i, j-1]:
                elem_new = (element[0], element[1]+".")
                minhp.add(elem_new)
            
            # Case 3
            for t in range(i, j):
                if match(s[t],  s[j]):
                        
                    # deal with number of structures
                    count_left = counts[i, t-1] if i<(t-1) else 1
                    count_right = counts[t+1, j-1] if (t+1)<(j-1) else 1
                    counts[i, j] += count_left*count_right
                    
                    # deal with best k
                    # structure on the left of t
                    seq_lefts = []
                    if i<=(t-1):
                        for seq_left in opts_k[i, t-1]:
                            seq_lefts.append((seq_left[0], seq_left[1]+"("))
                    else:
                        seq_lefts = [(0, "(")]
                    
                    # structure on the left of t
                    seq_rights = []
                    if (t+1)<=(j-1):
                        for seq_right in opts_k[t+1, j-1]:
                            seq_rights.append((seq_right[0], seq_right[1]+")"))
                    else:
                        seq_rights = [(0, ")")]
                    
                    # Cartesian product of left structures and right structures
                    for seq_left in seq_lefts:
                        for seq_right in seq_rights:
                            minhp.add(  (seq_left[0]+seq_right[0]+1, seq_left[1]+seq_right[1]) )
                    
            opts_k[i, j] = sorted(minhp.pool, reverse=True)
    return opts_k[0, -1][0], int(counts[0, -1]), opts_k[0, -1]