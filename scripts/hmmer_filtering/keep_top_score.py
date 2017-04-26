#! /usr/bin/env python
""" A script to keep only the best scoring hit for each gene. """

import argparse
import sys

def main(args):
    saved_rows = {}
    with open(args.infile) as ifh:
        for line in ifh.readlines():
            line_fields = line.split('\t')
            if line_fields[0] in saved_rows:
                # We know this line has a lower score, but does it have a lower e-value?
                if float(line_fields[4]) < saved_rows[line_fields[0]][0]:
                    sys.stderr.write(("E-value for gene {} inconsistent "
                    "with score, detected for profile {} "
                    "Evalues {} and {} \n").format(line_fields[0], line_fields[2], float(line_fields[4]), saved_rows[line_fields[0]][0]))
            else:
                saved_rows[line_fields[0]] = (float(line_fields[4]), line)
    for gene_id, eval_row_tuple in saved_rows.iteritems():
        sys.stdout.write(eval_row_tuple[1])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("infile", help=("Hmmsearch table input, "
            "where column 1 is gene id, column five is the e-value "
            "and column six is the score. Input file is sorted on score."))
    args = parser.parse_args()
    main(args)
