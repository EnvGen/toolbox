#!/usr/bin/env python
"""filter_bam_on_references.py

Filter a sam stream keeping the intended references and reads that either
map to the references or have a mate that map to the references. 

usage:

    samtools view -h input.bam | python --keep_all_mates filter_bam_on_references.py references.txt | samtools view -bS > output.bam
"""

import argparse
import sys
import pysam
import re

def is_mapped_to(read, references):
    read_fields = read.split('\t')
    if read_fields[2] in references:
        # If bit 4 is set, the read is unmapped
        return int(read_fields[1])&(1<<4) == 0

def read_reference_list(reference_file):
    references = []
    with open(reference_file, 'r') as ifh:
        for line in ifh:
            references.append(line.strip())
    return references


def main(args):
    header_re = re.compile("^@[A-Z][A-Z]")
    reference_header_re = re.compile("@SQ\tSN:([ -~]+)")
    references = read_reference_list(args.references)
    for line in sys.stdin:
        line = line.strip()
        # Header, output the ones we'd like to keep
        if header_re.match(line):
            m = reference_header_re.search(line)
            if m:
                if m.group(1) in references:
                    print(line)
            else:
                print(line)
        else:
            read1 = line
            break

    for i, line in enumerate(sys.stdin):
        line = line.strip()
        # Read two reads at a time
        if i % 2:
            read1 = line
        else:
            read2 = line
            both_printed=False
            if is_mapped_to(read1, references):
                if args.keep_all_mates:
                    print(read1, read2)
                    both_printed=True
                else:
                    print(read1)

            if not both_printed and is_mapped_to(read2, references):
                if args.keep_all_mates:
                    print(read1, read2)
                else:
                    print(read2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("references", help="References to keep in a file, one reference per line")
    parser.add_argument("--keep_all_mates", help="Keep all mates for mapping reads.")
    args = parser.parse_args()

    main(args)
