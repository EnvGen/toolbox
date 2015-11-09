#!/usr/bin/env python
"""filter_fasta_on_length.py 

Filter a fasta file so that the resulting sequences are all longer 
or equal to the length threshold parameter. Only accepts stdin."""

import argparse
import sys
from Bio import SeqIO

def main(args):
    seqs = []
    for i, seq in enumerate(SeqIO.parse(sys.stdin, "fasta")):
        if len(seq) >= args.length_threshold:
            seqs.append(seq)
        if i % 1000 == 0 and len(seqs):
            SeqIO.write(seqs, sys.stdout, 'fasta')
            seqs = []
    
    if len(seqs):
        SeqIO.write(seqs, sys.stdout, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-l", "--length_threshold", type=int, help="Length trheshold to filter on.")

    args = parser.parse_args()

    main(args)
