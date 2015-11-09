#!/usr/bin/env python
"""filter_fasta_on_length.py 

Filter a fasta file so that the resulting sequences are all longer 
or equal to the length threshold parameter."""

import argparse
import sys
from Bio import SeqIO

def main(args):
    seqs = []
    for i, seq in enumerate(SeqIO.parse(args.fasta_file, "fasta")):
        if len(seq) >= args.length_threshold:
            seqs.append(seq)
    if args.output_file:
        with open(args.output_file, 'w') as ofh:
            SeqIO.write(seqs, ofh, 'fasta')
    else:
        SeqIO.write(seqs, sys.stdout, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_file", help="Input Fasta file.")
    parser.add_argument("-l", "--length_threshold", type=int, help="Length trheshold to filter on.")
    parser.add_argument('-o', '--output_file',
        help=("Optional output file where sequences will be printed." 
            " Otherwise use stdout."))

    args = parser.parse_args()

    main(args)
