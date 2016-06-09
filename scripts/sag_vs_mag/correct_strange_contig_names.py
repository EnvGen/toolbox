#!/usr/bin/env python
"""A script to create plots for coverage for MAG contigs in SAG reads"""
import argparse
import pandas as pd
import sys
import gzip

def read_fasta_headers(fasta_file, compression=None):
    fasta_headers = []
    if compression == 'gzip':
        with gzip.open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    line = line.strip()
                    fasta_headers.append(line[1:].split(' ')[0])
    else:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    line = line.strip()
                    fasta_headers.append(line[1:].split(' ')[0])
    return fasta_headers

def main(args):
    regular_to_strange = {}
    to_be_corrected_df = pd.read_table(args.to_be_corrected, names=['contig_id'])
    fasta_headers = read_fasta_headers(args.correct_fasta_file, compression='gzip')

    wrong_to_correct = dict((fasta_header.split('_')[0], fasta_header) for fasta_header in fasta_headers)

    with open(args.output_file, 'w') as ofh:
        for wrong_contig_id in to_be_corrected_df['contig_id'].values:
            if wrong_contig_id in wrong_to_correct:
                ofh.write(wrong_to_correct[wrong_contig_id] + '\n')
            else:
                sys.stderr.write("Did not find corrected id for {}\n".format(wrong_contig_id))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('to_be_corrected', help='Single column file with contig ids')
    parser.add_argument('correct_fasta_file', help='A fasta file with correct contig ids')
    parser.add_argument('output_file', help='The new list of corrected contig ids')
    args = parser.parse_args()
    main(args)
