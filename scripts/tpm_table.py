#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
"""A script to calculate TPM values for contigs or genes based on coverage files

TPM values are defined as in Wagner et al (Theory in Biosciences) 2012. 

Except rg x rl / flg is replaced with average coverage.

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

        T = read_length * count_df['count'].divide(gene_lengths['gene_length']).sum()

        # Rpkm calculation based on average coverage
        tpm = ((1e6*read_length)/float(T))*(count_df['count'].divide(gene_lengths['gene_length']))
        
        df[sample_name] = tpm

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
