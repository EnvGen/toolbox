#!/usr/bin/env python
"""filter_fasta_on_length.py 

Filter a fasta file so that the resulting sequences are all longer 
or equal to the length threshold parameter. Only accepts stdin."""

import argparse
import sys
from Bio import SeqIO

def main(args):
    seqs = []
    if args.input_fasta is not None:
        input_handle = open(args.input_fasta, 'r')
    else:
        input_handle = sys.stdin

    for i, seq in enumerate(SeqIO.parse(input_handle, "fasta")):
        if (len(seq) >= args.length_threshold and not args.shorter) or (len(seq) < args.length_threshold and args.shorter):
            seqs.append(seq)
        if i % 1000 == 0 and len(seqs):
            SeqIO.write(seqs, sys.stdout, 'fasta')
            seqs = []
    
    if len(seqs):
        SeqIO.write(seqs, sys.stdout, 'fasta')

    input_handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_fasta", help="Input file")
    parser.add_argument("-l", "--length_threshold", type=int, help="Length threshold to filter on.")
    parser.add_argument("--shorter", action="store_true", help="Output contigs shorter than threshold, default is to output longer or equal sequences.")
    args = parser.parse_args()

    main(args)
