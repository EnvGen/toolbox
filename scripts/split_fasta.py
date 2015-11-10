#!/usr/bin/env python
"""A script to split a fasta file into multiple smaller files.

Output files will be <original>.xx.fasta
"""
from __future__ import print_function
import argparse
import sys
import os
import re

def header(row):
    regex = re.compile("^>")
    return regex.findall(row)

def next_output_file(input_file, prefix, i):
    if input_file is not None:
        prefix = ".".join(os.path.basename(input_file).split('.')[0:])
    
    new_filename = "{0}.{1}.fasta".format(prefix, i)
    print(new_filename)
    return open(new_filename, 'w')


def main(args):
    if args.input_fasta is not None:
        input_handle = open(args.input_fasta, 'r')
    else:
        input_handle = sys.stdin

    n = 0
    i = 0
    outputfh = next_output_file(args.input_fasta, args.prefix, i) 
    for row in input_handle:
        if header(row):
            n+=1

        if n > args.n_seqs:
            i += 1
            n = 1
            outputfh.close()
            outputfh = next_output_file(args.input_fasta, args.prefix, i)
        
        outputfh.write(row)

    outputfh.close()
    input_handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input_fasta', help="Input read fasta 1")
    parser.add_argument('--prefix', help="output prefix")
    parser.add_argument('n_seqs', type=int, help="Number of sequences per file, the last file will contain slightly less")
    args = parser.parse_args()
    main(args)
