#!/usr/bin/env python
"""A script to fetch all dbxref annotations of a certain type from a gff file"""

import argparse
import re

TYPE_REGS = { "PFAM": 'Dbxref=pfam([0-9]+)',
              "COG":  'Dbxref=COG([0-9]+)',
              "TIGR": 'Dbxref=TIGR([0-9]+)'}

TYPE_REGS_WITH_EVALUE = { "PFAM": 'Dbxref=pfam([0-9]+),evalue:([e\-\d\.]*),',
                          "COG": 'Dbxref=COG([0-9]+),evalue:([e\-\d\.]*),',
                          "TIGR": 'Dbxref=TIGR([0-9]+),evalue:([e\-\d\.]*),'}

def main(args):
    if args.with_evalue:
        reg_s = TYPE_REGS_WITH_EVALUE[args.type]
    else:
        reg_s = TYPE_REGS[args.type]

    reg = re.compile(reg_s)
    gene_reg = re.compile('ID=([^;]+)')

    skip_rest = False
    matches = []
    with open(args.infile, 'r') as fh:
        for line in fh:
            if skip_rest:
                break
            if line.startswith('>'):
                skip_rest = True
            elif line.startswith('#'):
                continue
            else:
                gene_match = gene_reg.search(line)
                if gene_match is None:
                    print(line)
                    exit(-1)
                gene_id = gene_match.group(1)
                for match in reg.findall(line):
                    if args.with_evalue:
                        print("\t".join([gene_id, str(args.type + match[0]), match[1]]))
                    else:
                        print("\t".join([gene_id, str(args.type + match)])) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", 
            help=("Gff file with annotations as dbxref"))
    parser.add_argument("type", choices=["PFAM", "TIGR", "COG", "EC"],
            help=("The annotation type to fetch"))
    parser.add_argument("--with_evalue", action='store_true')
    args = parser.parse_args()
    main(args)
