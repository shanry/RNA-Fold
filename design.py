import argparse
import os
import time

from nupack import *
import pandas as pd

my_model = Model(material='RNA', celsius=37)

def design_batch(structs_list):
    defect_list = []
    for i, struct in enumerate(structs_list):
        start = time.time()
        seq_des = des(structure=struct, model=my_model)
        end = time.time()
        df = defect(strands=seq_des, structure=struct, model=my_model)
        defect_list.append((len(struct), df, end-start))
        print(f'{i}, length: {len(struct)}, defect: {df:.4f}, time: {end-start:.2f}')
    return defect_list


def main(path):
    with open(path) as f:
        structs_list = [line.strip() for line in f.readlines()]
        results = design_batch(structs_list)
    df = pd.DataFrame(results, columns = ['length', 'defect', 'time'])
    df.to_csv(path+".csv")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, default="data/9families_free/tmRNA_ref")
    args = parser.parse_args()
    print(f"args: {args}")
    main(args.path)
    