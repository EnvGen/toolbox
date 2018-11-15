#!/usr/bin/env python
"""
Based on the busco results, approves bins according to the leves of contamination and completeness.

Copies approved bins to output directory.

@author: alneberg
"""
from __future__ import print_function
import sys
import os
import argparse
import pandas as pd
from shutil import copyfile

def main(args):
    # Read in the busco table
    df = pd.read_table(args.busco_stats, index_col=0, sep=',')
    df['completeness_float'] = df['completeness'].apply(lambda x: float(x.replace('%','')))
    df['duplicated_float'] = df['duplicated'].apply(lambda x: float(x.replace('%','')))
    # extract the ids for all rows that meet the requirements
    filtered_df = df[(df['completeness_float'] >= args.min_completeness) & (df['duplicated_float'] <= args.max_contamination)] 
    
    if len(filtered_df) < 1:
        sys.stderr.write("\nApproved 0 bins\n\n")
        return None

    approved_bins = list(filtered_df.index)

    # copy the approved bins to the new output directory
    for approved_bin_int in approved_bins:
        approved_bin = str(approved_bin_int)
        bin_source = os.path.join(args.bin_directory, approved_bin)
        bin_source += '.' + args.extension
        bin_destination = os.path.join(args.output_directory)
        bin_destination += '/' + os.path.basename(bin_source)
        
        sys.stderr.write("Copying approved bin {} from {} to {}\n".format(approved_bin, bin_source, bin_destination))
        copyfile(bin_source, bin_destination)

    sys.stderr.write("\nApproved {} bins\n\n".format(len(approved_bins)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bin_directory", help=("Input fasta files should be within directory."))
    parser.add_argument("busco_stats", help="Busco stats in csv format with header and bin id as index/first column")
    parser.add_argument("output_directory", help="Directory where to put approved bins")
    parser.add_argument("--min_completeness", default=30, type=float, help="default=30")
    parser.add_argument("--max_contamination", default=10, type=float, help="default=10")
    parser.add_argument("--extension", default='fa')
    args = parser.parse_args()

    main(args)
