#!/usr/bin/env python
from __future__ import print_function
"""A script to print cog categories given the cog ids per gene."""

import sys
import argparse
import pandas as pd

def main(args):
    category_df = pd.read_table(args.cog_category_list_file, sep=',', header=None, index_col=0, names=['cog_id', 'category', 'annotation_name'])
    gene_annotation_df = pd.read_table(args.gene_annotation_file, header=None, index_col=0, names=['gene_id', 'cog_id', 'evalue'])

    gene_annotation_df['cog_category'] = gene_annotation_df['cog_id'].apply(lambda x: category_df.ix[x]['category'])

    gene_annotation_df.to_csv(sys.stdout, sep='\t', header=None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_annotation_file",
            help=("TSV file with cog ids on "))
    parser.add_argument("cog_category_list_file",
            help=("TSV file with category for all COG ids"))
    args = parser.parse_args()
    main(args)
