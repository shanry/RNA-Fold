import argparse

from k_best import algo_0

def test(k):
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
    for rna in rnas:
        best_1, num_struct, best_k = algo_0(rna, k)
        print(rna)
        print(best_1)
        print(num_struct)
        print(best_k)
        print()
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", type=str, default="GCACG")
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--test", action='store_true')
    args = parser.parse_args()
    print(f"args: {args}")
    best_1, num_struct, best_k = algo_0(args.rna, args.k)
    print()
    print(args.rna)
    print(best_1)
    print(num_struct)
    print(best_k)
    
    if args.test:
        test(args.k)
    