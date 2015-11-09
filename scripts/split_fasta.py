#!/usr/bin/env python
"""A script to split a fasta file into multiple smaller files.

Output files will be <original>.xx.fasta
"""

import argparse
import sys
import os
import re

def header(row):
    regex = re.compile("^>")
    return regex.findall(row)

def next_output_file(input_file, i):
    orig_name = ".".join(os.path.basename(input_file).split('.')[0:])
    print "{0}.{1}.fasta".format(orig_name, i)
    return open("{0}.{1}.fasta".format(orig_name, i), 'w')


def main(args):
    with open(args.input_fasta, 'r') as input:
        n = 0
        i = 0
        outputfh = next_output_file(args.input_fasta, i) 
        for row in input:
            if header(row):
                n+=1

            if n > args.n_seqs:
                i += 1
                n = 1
                outputfh.close()
                outputfh = next_output_file(args.input_fasta, i)
            
            outputfh.write(row)

        outputfh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_fasta', help="Input read fasta 1")
    parser.add_argument('n_seqs', type=int, help="Number of sequences per file, the last file will contain slightly less")
    args = parser.parse_args()
    main(args)
