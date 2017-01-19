#! /usr/bin/env python
""" A script to keep hits that pass the strict criterion e-value < 1e-18 and coverage > 0.35
corresponding to the recommended criterion for bacteria  by dbCAN.

Given http://csbl.bmb.uga.edu/dbCAN/download/readme.txt
"""

import argparse
import sys

def main(args):
    domain_score = {}
    with open(args.domain_pos_and_score) as ifh:
        for line in ifh.readlines():
            line = line.strip()
            hmm_id, protein_id, e_value, aln_start, aln_end, score = line.split('\t')
            domain_score[(hmm_id, protein_id, e_value, aln_start, aln_end)] = score

    with open(args.infile) as ifh:
        for line in ifh.readlines():
            line = line.strip()
            line_fields = line.split('\t')

            hmm_id, protein_id, e_value, aln_start, aln_end = (line_fields[0], line_fields[2], line_fields[4], line_fields[5], line_fields[6])
            
            try:
                score = domain_score[(hmm_id, protein_id, e_value, aln_start, aln_end)]
            except KeyError:
                sys.stderr.write("Tuple: {}".format((hmm_id, protein_id, e_value, aln_start, aln_end)))
                raise

            if float(line_fields[4]) < float(1e-18):
                if float(line_fields[9]) > 0.35:
                    sys.stdout.write(line + '\t' + score + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("infile", help=("output from dbcan hmmscan-parser.sh script"))
    parser.add_argument("domain_pos_and_score", help="A file with the score for each domain alignment, contains positions to assert uniqueness")
    args = parser.parse_args()
    main(args)
