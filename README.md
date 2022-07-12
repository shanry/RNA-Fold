# RNA algorithms
Folding and partition of RNA
## RNA fold
K-best parsing for the secondary structure prediction of RNA. \
Implementation of algorithms from the paper [Better k-best Parsing](https://aclanthology.org/W05-1506.pdf).

### Algorithm 0 
python fold.py --rna AGGCAUCAAACCCUGCAUGGGAGCG --k 10 --algo 0 \
complexity: $\mathcal{O}(n^3\cdot k^2\log{k})$

### Algorithm 1 
python fold.py --rna AGGCAUCAAACCCUGCAUGGGAGCG --k 10 --algo 1 \
complexity: $\mathcal{O}(n^3\cdot k\log{k})$

### Algorithm 2
python fold.py --rna AGGCAUCAAACCCUGCAUGGGAGCG --k 10 --algo 2 \
complexity: $\mathcal{O}(n^3 + n^2\cdot k\log{k})$

### Algorithm 3
python fold.py --rna AGGCAUCAAACCCUGCAUGGGAGCG --k 10 --algo 3 \
complexity: $\mathcal{O}(n^3 + n\cdot k\log{k})$

### Get all secondary structures
python fold.py --rna UUGGACUUG --algo 4 # the output can be very long\
complexity: $\mathcal{O}(n^3)$

## Partition function (count version)
### Count of inside and outside (triple output)
python partition.py --inside forward  --outside forward  --sharpturn 0  --rna UUUGGCACUA \
python partition.py --inside forward  --outside backward --sharpturn 0 --rna UUUGGCACUA \
python partition.py --inside backward --outside forward  --sharpturn 0 --rna UUUGGCACUA \
python partition.py --inside backward --outside backward --sharpturn 0 --rna UUUGGCACUA