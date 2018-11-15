#!/usr/bin/env python
"""
Based on the checkm results, approves bins according to the leves of contamination and completeness.

Prints the corresponding stats to stdout.

@author: alneberg
"""
from __future__ import print_function
import sys
import os
import argparse
import pandas as pd
from shutil import copyfile

def main(args):
    # Read in the checkm table
    df = pd.read_table(args.checkm_stats, index_col=0)
    # extract the ids for all rows that meet the requirements
    filtered_df = df[(df['Completeness'] >= args.min_completeness) & (df['Contamination'] <= args.max_contamination)] 
    
    filtered_df.to_csv(sys.stdout, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkm_stats", help="Checkm qa stats in tab_table format")
    parser.add_argument("--min_completeness", default=85, type=float, help="default=85")
    parser.add_argument("--max_contamination", default=5, type=float, help="default=5")
    args = parser.parse_args()

    main(args)
