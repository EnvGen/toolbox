#!/usr/bin/env python
"""A script to create plots for coverage for MAG contigs in SAG reads"""
import argparse
import pandas as pd
import sys

def read_cov_file(cov_file):
    avg_cov_d = {}
    df = pd.read_table(cov_file, header=None, names=['contig_id', 'depth', 'bases_with_depth', 'tot_bases', 'fraction_of_bases_with_depth'])
    for contig_id, contig_df in df.groupby('contig_id'):
        try:
            fraction_not_cov = contig_df[contig_df['depth'] == 0]['fraction_of_bases_with_depth'].values[0]
        except IndexError:
            fraction_not_cov = 0

        bases_covered = contig_df[contig_df['depth'] != 0]['bases_with_depth'].sum()
        
        avg_cov_d[contig_id] = {'bases_covered': bases_covered, 
                'fraction_not_cov': fraction_not_cov, 
                'avg_cov': (contig_df['depth'] * contig_df['fraction_of_bases_with_depth']).sum(), 
                'length': contig_df['tot_bases'].max()}
    avg_cov_d.pop('genome')
    return pd.DataFrame.from_dict(avg_cov_d, orient='index')


def main(args):
    assert(len(args.input_coverage) == len(args.mag_contig_list))
    for input_coverage, contig_list in zip(args.input_coverage, args.mag_contig_list):
        avg_cov_df = read_cov_file(input_coverage)
        bin_contigs = pd.read_table(contig_list, names=['contig_id'])
        bin_contigs_s = set(bin_contigs["contig_id"])
        avg_cov_df['in_bin'] = [contig in bin_contigs_s for contig in avg_cov_df.index]

        #avg_cov_df[avg_cov_df['in_bin']]
        print(input_coverage, contig_list)
        avg_cov_df[(~avg_cov_df['in_bin']) & (avg_cov_df['length'] > 1000)][['avg_cov', 'bases_covered', 'fraction_not_cov']].to_csv(sys.stdout, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_coverage', nargs='*', help='Bedtools coverage files')
    parser.add_argument('--mag_contig_list', nargs='*', help='contig list with all contigs in MAG')
    args = parser.parse_args()
    main(args)
