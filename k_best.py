import numpy as np

import heapq

pairs = {'au', 'gc', 'gu'}

def match(x, y):
    return (x+y).lower() in pairs or (y+x).lower() in pairs


class MinHeap:
 
    def __init__(self, k):
        self.pool = []
        self.k = k
 
    def add(self, item):
 
        # if the min-heap's size is less than `k`, push directly
        if len(self.pool) < self.k:
            heapq.heappush(self.pool, item)
            return True
 
        # otherwise, compare to decide push or not
        elif self.pool[0][0] < item[0]:
            heapq.heappushpop(self.pool, item)
            return True
        return False
    
    def pop(self):
        return heapq.heappop(self.pool)
            

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


def algo_1(s, k=10):
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
            for element in opts_k[i, j-1][:k]:
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
                        for seq_left in opts_k[i, t-1][:k]:
                            seq_lefts.append((seq_left[0], seq_left[1]))
                    # else:
                    #     seq_lefts = [(0, "(")]
                    
                    # structure on the left of t
                    seq_rights = []
                    if (t+1)<=(j-1):
                        for seq_right in opts_k[t+1, j-1][:k]:
                            seq_rights.append((seq_right[0], seq_right[1]))
                    # else:
                    #     seq_rights = [(0, ")")]
                    
                    if not seq_lefts and not seq_rights:
                        minhp.add(  (1, '()' ) )
                    elif not seq_lefts and seq_rights:
                        for seq_right in seq_rights:
                            minhp.add(  (seq_right[0]+1, '('+seq_right[1]+')') )
                    elif seq_lefts and not seq_rights:
                        for seq_left in seq_lefts:
                            minhp.add(  (seq_left[0]+1, seq_left[1]+'()') )
                    else:
                        # item_start = (seq_lefts[0][0]+seq_rights[0][0]+1, seq_lefts[0][1]+'('+seq_rights[0][1]+')') 
                        bitset = set() # indicate the existence of row*k+col
                        mincand = MinHeap(k*k)
                        cand = (-seq_lefts[0][0]-seq_rights[0][0]-1, seq_lefts[0][1]+'('+seq_rights[0][1]+')', 0, 0)
                        mincand.add(cand)
                        bitset.add(0*k+0)
                        while len(mincand.pool)>0:
                            cand_best = mincand.pop()
                            # print('cand_best: ', cand_best)
                            item_best = (-cand_best[0], cand_best[1])
                            if minhp.add(item_best):
                                loc = (cand_best[-2], cand_best[-1])
                                # axis left
                                loc_left = loc[0]+1
                                loc_right = loc[1]
                                if loc_left < len(seq_lefts) and loc_left*k+loc_right not in bitset:
                                    cand = (-seq_lefts[loc_left][0]-seq_rights[loc_right][0]-1, seq_lefts[loc_left][1]+'('+seq_rights[loc_right][1]+')', loc_left, loc_right)
                                    mincand.add(cand)
                                    bitset.add(loc_left*k+loc_right)
                                # axis right
                                loc_left = loc[0]
                                loc_right = loc[1]+1
                                if loc_right < len(seq_rights) and loc_left*k+loc_right not in bitset:
                                    cand = (-seq_lefts[loc_left][0]-seq_rights[loc_right][0]-1, seq_lefts[loc_left][1]+'('+seq_rights[loc_right][1]+')', loc_left, loc_right)
                                    mincand.add(cand)
                                    bitset.add(loc_left*k+loc_right)
                            else:
                                break
                            
                        # Cartesian product of left structures and right structures
                        # for seq_left in seq_lefts:
                        #     for seq_right in seq_rights:
                        #         minhp.add(  (seq_left[0]+seq_right[0]+1, seq_left[1]+'('+seq_right[1]+')') )
            opts_k[i, j] = sorted(minhp.pool, reverse=True)
    return opts_k[0, -1][0], int(counts[0, -1]), opts_k[0, -1]