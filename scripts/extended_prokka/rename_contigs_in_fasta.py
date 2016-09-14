#!/usr/bin/env python
"""rename_contigs_in_fasta.py 

Append string to contig id for each contig in input fasta file.
Prints to stdout.
"""

import argparse
import sys
from Bio import SeqIO

def main(args):
    seqs = []
    for seq in SeqIO.parse(args.fasta_file, "fasta"):
        new_id = seq.id + "_" + args.string_to_append
        seq.description = seq.description.replace(seq.id, new_id)
        seq.id = new_id
        SeqIO.write([seq], sys.stdout, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_file", help="Input Fasta file.")
    parser.add_argument("-s", "--string_to_append", help="String to append to each contig")

    args = parser.parse_args()

    main(args)
