#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
"""A script to calculate rpkm values for contigs or genes based on coverage files

Output:
    Tab separated values: gene id, average coverage and gene length, printed to stdout.

"""
import sys
import argparse
import pandas as pd

def main(args):
    sample_info = pd.read_table(args.sample_info, header=None, index_col=0)
    df = pd.DataFrame()
    gene_lengths = pd.read_table(args.gene_lengths, header=None, index_col=0, names=['gene_id', 'gene_length'])

    for fn, sample_name in zip(args.coverage_files, args.sample_names):
        count_df = pd.read_table(fn, index_col=0, header=None, names=['gene_id', 'count'])
        nr_reads_sample, read_length = sample_info.ix[sample_name]

        # Rpkm calculation based on average coverage
        rpkm = (count_df['count'].divide(gene_lengths['gene_length'])).divide(float(nr_reads_sample)) * 1e9
        
        df[sample_name] = rpkm
    df.to_csv(sys.stdout, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sample_names', nargs='*', 
            help=("Sample names, in the same order as coverage_files"))
    parser.add_argument('--coverage_files', nargs='*', 
            help=("Coverage files with tab separated values: "
                "sequence id, average coverage, sequence length"))
    parser.add_argument('--sample_info', 
            help=("Tab separated values 'sample_id', 'nr_reads', 'avg_read_length'. "
                "all values in sample_names need to be present as sample_id values"))
    parser.add_argument('--gene_lengths',
            help=("Gene lengths in a tsv file"))
    args = parser.parse_args()
    main(args)
