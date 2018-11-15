#!/usr/bin/env python
"""
Reads an augustus output file and prints the amino_acid fasta file to stdout

@author: alneberg
"""
import sys
import os
import argparse

def to_fasta(s, gene_id, contig_id):
    print('>{}_{}'.format(contig_id, gene_id))
    print(s)
    

def main(args):
    with open(args.augustus_output) as ifh:
        protein_str = None
        for line in ifh: 
            line = line.strip()
            # Check contig_id
            if not line.startswith('#'):
                contig_id = line.split('\t')[0]
            # Check if protein starts
            elif line.startswith('# protein sequence'):
                line = line.replace('# protein sequence = [', '')
                protein_str = line.replace(']','')
                protein_str = line
            elif protein_str:
                # If protein ends
                # Parse lines and output to stdout
                if line.startswith('# end gene'):
                    line = line.replace('# end gene ', '')
                    gene_id = line
                    to_fasta(protein_str, gene_id, contig_id)
                    protein_str = None
                else:
                    # add to protein lines
                    line = line.replace('# ', '')
                    line = line.replace(']', '')
                    protein_str += line

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("augustus_output", help=("Standard output format of augustus."))
    args = parser.parse_args()

    main(args)
