#!/usr/bin/env python
"""A script to create plots for coverage for MAG contigs in SAG reads"""
import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def read_cov_file(cov_file):
    avg_cov_d = {}
    df = pd.read_table(cov_file, header=None, names=['contig_id', 'depth', 'bases_with_depth', 'tot_bases', 'fraction_of_bases_with_depth'])
    for contig_id, contig_df in df.groupby('contig_id'):
        avg_cov_d[contig_id] = {'avg_cov': (contig_df['depth'] * contig_df['fraction_of_bases_with_depth']).sum(), 'length': contig_df['tot_bases'].max()}
    avg_cov_d.pop('genome')
    return pd.DataFrame.from_dict(avg_cov_d, orient='index')


def main(args):
    assert(len(args.input_coverage) == len(args.mag_contig_list))
    ax = plt.subplot()
    for input_coverage, contig_list in zip(args.input_coverage, args.mag_contig_list):
        avg_cov_df = read_cov_file(input_coverage)
        bin_contigs = pd.read_table(contig_list, names=['contig_id'])
        bin_contigs_s = set(bin_contigs["contig_id"])
        avg_cov_df['in_bin'] = [contig in bin_contigs_s for contig in avg_cov_df.index]

        avg_cov_df[avg_cov_df['in_bin']].plot.scatter(x='avg_cov', y='length', color='green', ax=ax, alpha=args.alpha)
        avg_cov_df[(~avg_cov_df['in_bin'])].plot.scatter(x='avg_cov', y='length', color='red', ax=ax, alpha=args.alpha)
        ax.set_xlim((0.001, 100000))
        ax.set_ylim((100, 1000000))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel("Coverage Depth for SAG reads")
    plt.savefig(args.output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_coverage', nargs='*', help='Bedtools coverage files')
    parser.add_argument('--mag_contig_list', nargs='*', help='contig list with all contigs in MAG')
    parser.add_argument('--output_file', default='coverage_plot.pdf')
    parser.add_argument('--alpha', default=0.7, help='The value for alpha transparancy in the scatter plot')
    args = parser.parse_args()
    main(args)
