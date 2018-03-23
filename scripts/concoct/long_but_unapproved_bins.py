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
from shutil import copyfile


def main(args):
    # Read bin sizes file
    bin_sizes_df = pd.read_table(args.bin_sizes, sep=',', index_col=0)

    if os.stat(args.approved_bins).st_size == 0:
        approved_bins_df = None
    else:
        approved_bins_df = pd.read_table(args.approved_bins, header=None, index_col=0) 

    for bin_nr, row in bin_sizes_df.iterrows():
        if float(row['nr_bases']) < 1000000:
            continue

        # create a symlink unless bin is approved
        if (approved_bins_df is None) or (bin_nr not in approved_bins_df.index):
            source_bin = os.path.join(args.input_dir, "{}.fa".format(bin_nr))
            destination_bin = os.path.join(args.output_dir, "{}.fa".format(bin_nr))
            copyfile(source_bin, destination_bin)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir")
    parser.add_argument("bin_sizes")
    parser.add_argument("approved_bins")
    parser.add_argument("output_dir")

    args = parser.parse_args()
    main(args)
