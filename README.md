# RNA-Fold
K-best parsing for the secondary structure prediction of RNA. \
Implementation of algorithms from the paper [Better k-best Parsing](https://aclanthology.org/W05-1506.pdf).

## Algorithm 0 
python fold.py --rna AGGCAUCAAACCCUGCAUGGGAGCG --k 10 --algo 0 \
complexity: $\mathcal{O}(n^3\cdot k^2\log{k})$

## Algorithm 1 
python fold.py --rna AGGCAUCAAACCCUGCAUGGGAGCG --k 10 --algo 1 \
complexity: $\mathcal{O}(n^3\cdot k\log{k})$

## Algorithm 2
python fold.py --rna AGGCAUCAAACCCUGCAUGGGAGCG --k 10 --algo 2 \
complexity: $\mathcal{O}(n^3 + n^2\cdot k\log{k})$

## Algorithm 3
python fold.py --rna AGGCAUCAAACCCUGCAUGGGAGCG --k 10 --algo 3 \
complexity: $\mathcal{O}(n^3 + n\cdot k\log{k})$