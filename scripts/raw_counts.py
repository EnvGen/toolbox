#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
"""A script to combine raw counts for contigs or genes

"""
import sys
import argparse
import pandas as pd

def main(args):
    gene_lengths = pd.read_table(args.gene_lengths, header=None, index_col=0, names=['gene_id', 'gene_length'])

    df = None
    for fn, sample_name in zip(args.coverage_files, args.sample_names):
        count_df = pd.read_table(fn, index_col=0, header=None, 
                names=['gene_id', sample_name], compression=args.input_compression)
        if df is None:
            df = count_df
        else:
            df[sample_name] = count_df[sample_name]

    df['gene_length'] = gene_lengths
    df.to_csv(sys.stdout, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sample_names', nargs='*', 
            help=("Sample names, in the same order as coverage_files"))
    parser.add_argument('--coverage_files', nargs='*', 
            help=("Coverage files with tab separated values: "
                "sequence id, average coverage, sequence length"))
    parser.add_argument('--gene_lengths',
            help=("Gene lengths in a tsv file"))
    parser.add_argument("--input_compression", default=None, choices=[None, 'gzip'], 
            help="Compression type for input coverage files. Default=None, use 'gzip', for gzipped files.")
    args = parser.parse_args()
    main(args)
