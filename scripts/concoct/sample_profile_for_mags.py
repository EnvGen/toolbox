#!/usr/bin/env python
"""
With contigs cutup with cut_up_fasta.py as input, sees to that the consequtive
parts of the original contigs are merged.

prints result to stdout.

@author: alneberg
"""
import sys
import os
import argparse
import pandas as pd

def original_contig_name(s):
    """Transform s to the original contig name"""
    n = s.split(".")[-1]
    try:
        int(n)
    except:
        return s
    # Only small integers are likely to be 
    # indicating a cutup part.
    if int(n) < 1000:
        return ".".join(s.split(".")[:-1])
    else:
        # A large n indicates that the integer
        # was part of the original contig
        return s


def main(args):

    # Check if any approved bins exists:
    if os.stat(args.list_all_approved_bins).st_size == 0:
        return

    print("# Find all approved bins")
    approved_bins = pd.read_table(args.list_all_approved_bins, header=None, index_col=0)

    print("# Open clustering file (non-cutup but")
    clustering_df = pd.read_table(args.clustering_nocutup, header=None, index_col=0, names=['contig_id', 'cluster_id'], sep=',')

    print("# Read input table")
    input_table = pd.read_table(args.concoct_inputtable, index_col=0)
    original_columns = input_table.columns

    print("# Read kallisto quant result for the contig length")
    kallisto_df = pd.read_table(args.kallisto_quant, index_col=0)

    print("# Create a new column with non-cutup contig id")
    input_table['target_id'] = input_table.index
    input_table['original_contig_id'] = input_table['target_id'].apply(original_contig_name)

    print("# Subset clustering file to get all approved contigs.")
    approved_contigs = clustering_df[clustering_df['cluster_id'].isin(approved_bins.index)]


    print("# Subset input table to all contigs included in approved bins")
    approved_table = input_table[input_table['original_contig_id'].isin(approved_contigs.index)]

    print("# Add length for the relevant contigs")
    approved_table['length'] = kallisto_df.loc[approved_table.index]['length']

    def cluster_id_for_approved(or_id):
        return approved_contigs.loc[or_id]['cluster_id']

    print("# Assign cluster id to input table")
    approved_table['cluster_id'] = approved_table['original_contig_id'].apply(lambda x: cluster_id_for_approved(x))


    results_d = {}

    # Group by cluster id
    for cluster, cluster_df in approved_table.groupby('cluster_id'):
        total_length = cluster_df['length'].sum() 
        cluster_name = "{}_bin{}".format(args.sample_name, cluster)
        cluster_df[original_columns] = cluster_df[original_columns].mul(cluster_df['length'], axis=0)
        result_s = cluster_df[original_columns].sum() / float(total_length)
        result_s['length'] = int(total_length)
        results_d[cluster_name] = result_s

    pd.DataFrame.from_dict(results_d, orient='index').to_csv(args.output_file, sep='\t', index_label="MAG")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("list_all_approved_bins", help=("A file with all approved bins."))
    parser.add_argument("clustering_nocutup", help="Clustering file for all no-cutup contigs")
    parser.add_argument("concoct_inputtable", help="concoct inputtable for the cutup contigs")
    parser.add_argument("kallisto_quant", help="quant result for one sample, to get contig length")
    parser.add_argument("sample_name", help="Sample name")
    parser.add_argument("output_file", help="output file")
    args = parser.parse_args()

    main(args)
