#! /usr/bin/env python
""" A script to keep hits that pass the strict criterion e-value < 1e-18 and coverage > 0.35
corresponding to the recommended criterion for bacteria  by dbCAN.

Given http://csbl.bmb.uga.edu/dbCAN/download/readme.txt
"""

import argparse
import sys

def main(args):
    saved_rows = {}
    with open(args.infile) as ifh:
        for line in ifh.readlines():
            line_fields = line.split('\t')
            if float(line_fields[4]) < float(1e-18):
                if float(line_fields[9]) > 0.35:
                    sys.stdout.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument("infile", help=("output from dbcan hmmscan-parser.sh script"))
    args = parser.parse_args()
    main(args)
