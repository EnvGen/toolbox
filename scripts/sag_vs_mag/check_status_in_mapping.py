#!/usr/bin/env python
usage="""Checks the status for a subset of reads within a bam file.

Three categories:
    Matching
    Soft clipped
    Not mapping
"""

import argparse
import sys
import re

soft_clip_re = re.compile('([0-9]+)S')
matched_re = re.compile('([0-9]+)M')

nr_mismatched_re = re.compile('XM:i:([0-9]+)')

def parse_cigar(cigar):
    soft_clipped_matches = soft_clip_re.findall(cigar)
    matched_matches = matched_re.findall(cigar)    
    soft_clipped_bases = sum([int(x) for x in soft_clipped_matches])
    matched_bases = sum([int(x) for x in matched_matches])

    return soft_clipped_bases, matched_bases

def main(args):
    all_reads = set()
    for line in open(args.subset_reads, 'r'):
        read = line.strip()
        all_reads.add(read)

    not_mapped = set()
    soft_clipped_mapped = set()
    other = set()
    good_mapped = set()

    # Read the input sam files
    for line in sys.stdin:
        line = line.strip()

        # Check if in header
        if line[0] != '@':
            # Not in header

            line_split = line.split('\t')
            if (int(line_split[1]) / 128) == 1:
                read_id = line_split[0] + "_2"
            else:
                read_id = line_split[0] + "_1"

            # Ignore reads which are not in the subset
            if read_id in all_reads:

                no_mapping_bit = (((((int(line_split[1]) % 128) % 64) % 32) % 16) % 8) / 4
                if no_mapping_bit == 1:
                    not_mapped.add(read_id)
                else:
                    cigar = line_split[5]

                    soft_clipped, matching = parse_cigar(cigar)
                    read_length = len(line_split[9])
                    mismatches = int(nr_mismatched_re.findall(line)[0])
                    
                    if soft_clipped > 20:
                        soft_clipped_mapped.add(read_id)
                    elif (matching > read_length*0.90):
                        good_mapped.add(read_id)
                    else:
                        other.add(read_id)

    print("Not Mapped\tSoft Clipped\tOther\tGood Mapped")
    print("{}\t{}\t{}\t{}".format(len(not_mapped), len(soft_clipped_mapped), len(other), len(good_mapped)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('subset_reads')
    args = parser.parse_args()

    main(args)
