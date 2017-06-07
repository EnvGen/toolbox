#!/usr/bin/env python
"""fasta_lengths.py 

prints the id and length of each sequence in input fasta.

If no input file is given, reads from stdin.
"""

import argparse
import sys
from Bio import SeqIO

def main(args):
    if args.input_fasta is not None:
        input_handle = open(args.input_fasta, 'r')
    else:
        input_handle = sys.stdin

    for i, seq in enumerate(SeqIO.parse(input_handle, "fasta")):
        sys.stdout.write("{}{}{}\n".format(seq.id, args.separator, len(seq)))
    
    input_handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_fasta", help="Input file")
    parser.add_argument("-s", "--separator", default='\t', type=str,\
            help="Character(s) to use as a delimiter. Default=\\t")

    args = parser.parse_args()

    main(args)
