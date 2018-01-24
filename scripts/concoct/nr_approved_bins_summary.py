#!/usr/bin/env python
"""
Based on all checkm results, creates a table containing the nr of approved genomes 
for all binning runs it can find within binning/*.

@author: alneberg
"""
from __future__ import print_function
import sys
import os
import argparse
import pandas as pd
import glob

def find_checkm_dirs():
    all_runs = {}
    for path in glob.glob("binning/*/*/output_*/*/checkm_output/stats.tsv"):
        run_d = {}
        path_parts = path.split('/')
        run_d["binner"] = path_parts[1]
        run_d["sample"] = path_parts[2]
        run_d["quant"] = "_".join(path_parts[3].split('_')[1:])
        run_d["run_params"] = path_parts[4]
        run_d['SpeedUp'] = 'SpeedUp_Mp' in path_parts[4]
        run_d['standardize'] = 'standardize' in path_parts[4]
        all_runs[path] = run_d
    return all_runs

def main(args):
    all_runs = find_checkm_dirs()

    for path, run_d in all_runs.items():
        # Read in the checkm table
        df = pd.read_table(path, index_col=0)
        # extract the ids for all rows that meet the requirements
        nr_approved = len(df[(df['Completeness'] >= args.min_completeness) & (df['Contamination'] <= args.max_contamination)])
        run_d['nr_approved'] = nr_approved

    result_df = pd.DataFrame.from_dict(all_runs, orient='index')
    result_df.to_csv(sys.stdout, sep='\t', columns=['binner', 'sample', 'quant', 'run_params', 'SpeedUp', 'standardize', 'nr_approved'], index_label='path')
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min_completeness", default=75, type=float, help="default=75")
    parser.add_argument("--max_contamination", default=5, type=float, help="default=5")
    args = parser.parse_args()
    main(args)
