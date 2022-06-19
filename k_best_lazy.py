import numpy as np

import heapq

pairs = {'au', 'gc', 'gu'}


def match(x, y):
    return (x+y).lower() in pairs or (y+x).lower() in pairs


def hash_grid(row, col, n):
    return n*row + col


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
        elif self.pool[0].score < item.score:
            heapq.heappushpop(self.pool, item)
            return True
        return False
    
    def pop(self):
        return heapq.heappop(self.pool)
    
    
class Cand:
    
    def __init__(self, start, end, left, right, score, label, idx):
        self.start = start
        self.end = end
        self.left = left
        self.right = right
        self.score = score
        self.label = label
    
    def __lt__(self, other):
        return self.score < other.score
    
    def __gt__(self, other):
        return self.score > other.score
    
    def __le__(self, other):
        return self.score <= other.score
    
    def __ge__(self, other):
        return self.score >= other.score
    
    def __eq__(self, other):
        return self.score == other.score
    
    
class HyperArc:
    
    def __init__(self, start, end, left, right, split, opt, k):
        self.start = start
        self.end = end
        self.left = left
        self.right = right
        self.split = split
        self.k = k
        self.opt = opt
        # self.seq_lefts = seq_lefts
        # self.seq_rights = seq_rights
        self.bitset = set()  # indicate the existence of row*k+col
        self.mincand = MinHeap(k*3)
        
    def initilize(self):
        mincand = self.mincand
        # seq_lefts = self.seq_lefts
        # seq_rights = self.seq_rights
        bitset = self.bitset
        k = self.k
        if not seq_lefts and not seq_rights:
            loc_left = -1
            loc_right = -1
            cand =  (-1, '()', loc_left, loc_right )
        elif not seq_lefts and seq_rights:
            loc_left = -1
            loc_right = 0
            cand = (-seq_rights[0][0]-1, '('+seq_rights[0][1]+')', loc_left, loc_right)
        elif seq_lefts and not seq_rights:
            loc_left = 0
            loc_right = -1
            cand = (-seq_lefts[0][0]-1, seq_lefts[0][1]+'()', loc_left, loc_right) 
        else:
            # bitset = set() # indicate the existence of row*k+col
            # mincand = MinHeap(k*3)
            loc_left = 0
            loc_right = 0
            cand = (-seq_lefts[0][0]-seq_rights[0][0]-1, seq_lefts[0][1]+'('+seq_rights[0][1]+')',loc_left, loc_right)
        self.mincand.add(cand)
        self.bitset.add(loc_left*k+loc_right)
        
    def extract_best(self):
        mincand = self.mincand
        seq_lefts = self.seq_lefts
        seq_rights = self.seq_rights
        bitset = self.bitset
        k = self.k
        if len(mincand.pool)>0:
            cand_best = mincand.pop()
            # print('cand_best: ', cand_best)
            item_best = (-cand_best[0], cand_best[1])
            # if minhp.add(item_best):
            loc = (cand_best[-2], cand_best[-1])
            # axis left
            loc_left = loc[0]+1
            loc_right = loc[1]
            if loc_left < len(seq_lefts) and loc_left*k+loc_right not in bitset:
                if loc_right == -1:
                    cand = (-seq_lefts[loc_left][0]-1, seq_lefts[loc_left][1]+'()', loc_left, loc_right)
                else:
                    cand = (-seq_lefts[loc_left][0]-seq_rights[loc_right][0]-1, seq_lefts[loc_left][1]+'('+seq_rights[loc_right][1]+')', loc_left, loc_right)
                mincand.add(cand)
                bitset.add(loc_left*k+loc_right)
            # axis right
            loc_left = loc[0]
            loc_right = loc[1]+1
            if loc_right < len(seq_rights) and loc_left*k+loc_right not in bitset:
                if loc_left == -1:
                    cand = (-seq_rights[loc_right][0]-1, '('+seq_rights[loc_right][1]+')', loc_left, loc_right)
                else:
                    cand = (-seq_lefts[loc_left][0]-seq_rights[loc_right][0]-1, seq_lefts[loc_left][1]+'('+seq_rights[loc_right][1]+')', loc_left, loc_right)
                mincand.add(cand)
                bitset.add(loc_left*k+loc_right)
            # else:
                # break
            return item_best
        else:
            return None
        

class HyperArcList:
    
    def __init__(self, hyperarc_list, k):
        self.k = k
        self.hyperarc_list = hyperarc_list
        self.mincand = MinHeap(k*3+len(self.hyperarc_list))
        
    def initilize(self):
        mincand = self.mincand
        for i, harc in enumerate(self.hyperarc_list):
            bitset = harc.bitset
            k = harc.k
            cand_first = Cand(harc.start, harc.end, harc.left, harc.right, score=-harc.opt[0], label=harc.opt[1], idx=i)
            mincand.add(cand_first)
            bitset.add(hash_grid(harc.left, harc.right, k*2))
        
    def extract_best(self):
        mincand = self.mincand
        if mincand.pool:
            # pop out the best
            cand_best = mincand.pop()
            # item_best = (-cand_best.score, cand_best.label)
            # # append next
            # idx_harc = cand_best[-1]
            # hyperarc = self.hyperarc_list[idx_harc]
            # seq_lefts = hyperarc.seq_lefts
            # seq_rights = hyperarc.seq_rights
            # k = hyperarc.k
            # bitset = hyperarc.bitset
            # # two axis
            # loc = (cand_best[2], cand_best[3])
            # # axis left
            # loc_left = loc[0]+1
            # loc_right = loc[1]
            # if loc_left < len(seq_lefts) and loc_left*k+loc_right not in bitset:
            #     if loc_right == -1:
            #         cand = (-seq_lefts[loc_left][0]-1, seq_lefts[loc_left][1]+'()', loc_left, loc_right, idx_harc)
            #     else:
            #         cand = (-seq_lefts[loc_left][0]-seq_rights[loc_right][0]-1, seq_lefts[loc_left][1]+'('+seq_rights[loc_right][1]+')', loc_left, loc_right, idx_harc)
            #     mincand.add(cand)
            #     bitset.add(loc_left*k+loc_right)
            # # axis right
            # loc_left = loc[0]
            # loc_right = loc[1]+1
            # if loc_right < len(seq_rights) and loc_left*k+loc_right not in bitset:
            #     if loc_left == -1:
            #         cand = (-seq_rights[loc_right][0]-1, '('+seq_rights[loc_right][1]+')', loc_left, loc_right, idx_harc)
            #     else:
            #         cand = (-seq_lefts[loc_left][0]-seq_rights[loc_right][0]-1, seq_lefts[loc_left][1]+'('+seq_rights[loc_right][1]+')', loc_left, loc_right, idx_harc)
            #     mincand.add(cand)
            #     bitset.add(loc_left*k+loc_right)
            # else:
                # break
            return cand_best
        else:
            return None
            



def algo_2(s, k=10):
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
            harc_list = []
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

                    # structure on the right of t
                    seq_rights = []
                    if (t+1)<=(j-1):
                        for seq_right in opts_k[t+1, j-1][:k]:
                            seq_rights.append((seq_right[0], seq_right[1]))
                    
                    harc = HyperArc(split=t, k=k, seq_lefts=seq_lefts, seq_rights=seq_rights)
                    harc_list.append(harc)
            hyperarc_list = HyperArcList(harc_list, k)
            hyperarc_list.initilize()
            while hyperarc_list.mincand.pool:
                if not minhp.add(hyperarc_list.extract_best()):
                    break
            opts_k[i, j] = sorted(minhp.pool, reverse=True)
    return opts_k[0, -1][0], int(counts[0, -1]), opts_k[0, -1]


def algo_3(s, k=10):
    assert len(s) > 1, "the length of rna should be at least 2!"
    n = len(s)
    bests = {}
    
    harcs = {}
    counts = np.ones((n, n)) # number of structures
    opts_k = np.zeros((n, n) ,dtype=object)
    for row in range(n):
        opts_k[row, row] = (0, '.') # tuple: (num_pair, prediction)
    # for row in range(n):
    #     for col in range(n):
    #         opts_k[row, col] = []
    #         if row == col:
    #             opts_k[row, col].append((0, '.'))  # tuple: (num_pair, prediction)
    for l in range(1, n):
        for i in range(0, n-l):
            j = i + l
            key = hash_grid(i, j, n)
            if key not in bests:
                bests[key] = []
            # Case 1 doesn't exist here because l starts from 1
            # Case 2
            counts[i, j] = counts[i, j-1]    
            opt_2 = opts_k[i, j-1][0], opts_k[i, j-1][1]+"."
            opts_k[i, j] = opt_2
            harc = HyperArc(start=i, end=j, left=1, right=-1, split=j, opt=opt_2, k=k)
            harc_list = []
            harc_list.append(harc)
            # Case 3
            for t in range(i, j):
                if match(s[t],  s[j]):
                        
                    # deal with number of structures
                    count_left = counts[i, t-1] if i<(t-1) else 1
                    count_right = counts[t+1, j-1] if (t+1)<(j-1) else 1
                    counts[i, j] += count_left*count_right
                    
                    # deal with best k
                    # structure on the left of t
                    if i<=(t-1):
                        seq_left = opts_k[i, t-1]
                        left = 1
                    else:
                        seq_left = (0, "")
                        left = 0

                    # structure on the right of t
                    if (t+1)<=(j-1):
                        seq_right = opts_k[t+1, j-1]
                        right = 1
                    else:
                        seq_right = (0, "")
                        right = 0
                    opt_1 = seq_left[0] + seq_right[0] + 1, seq_left[1] + "(" + seq_right[1] + ")"
                    harc = HyperArc(start=i, end=j, left=left, right=right, split=t, opt=opt_1, k=k)
                    harc_list.append(harc)
                    if opt_1[0] > opts_k[i, j][0]:
                        opts_k[i, j] = opt_1
            hyperarc_list = HyperArcList(harc_list, k)
            hyperarc_list.initilize()
            if key not in harcs:
                harcs[key] = hyperarc_list
            if hyperarc_list.mincand.pool:
                bests[key].append(hyperarc_list.extract_best())
                    
    return opts_k[0, -1], int(counts[0, -1]), bests[hash_grid(0, n-1, n)]


def lazy_kth_best(start, end, n, kth, k, bests, harcs):
    key = hash_grid(start, end, n)
    while len(bests[key]) < kth:
        if bests[key]:
            best_last = bests[key][-1]
            lazy_next(best_last, k, bests, harcs)
        if harcs[key].mincand.pool:
            bests[key].append(harcs[key].extract_best())
        else:
            break
    return key


def lazy_next(best_last, k, bests, harcs):
    key = hash_grid(best_last.left, best_last.right)
    harc = harcs[key]
    if best_last.left > 0:
        kth = best_last.left+1
        lazy_kth_best(harc.start, harc.split-1, n, kth,  k, bests, harcs)
        if kth <= len(bests[hash_grid(harc.start, harc.split-1, n)]) and hash_grid(kth, best_last.right) not in harc.harc[best_last.idx]:
            bests[key].append(bests[hash_grid(harc.start, harc.split-1, n)][kth-1])
    return

