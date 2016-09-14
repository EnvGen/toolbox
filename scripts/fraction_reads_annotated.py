#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
"""A script to subset raw counts on annotated genes

"""
import sys
import argparse
import pandas as pd

def main(args):
    annotated_genes_df = pd.read_table(args.annotated_genes, header=None, index_col=0)

    count_df = pd.read_table(args.raw_counts, index_col=0)

    val_cols = list(count_df.columns)
    val_cols.remove('gene_length')
    subset_df = count_df.ix[list(annotated_genes_df.index)][val_cols]
    sys.stdout.write(str(subset_df.sum()))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--raw_counts',
            help=("File with tab separated values of raw counts "
                "per gene and a gene_length column"))
    parser.add_argument('--annotated_genes',
            help=("Annotated genes"))
    args = parser.parse_args()
    main(args)
