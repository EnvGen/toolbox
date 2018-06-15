#!/usr/bin/env python
usage="""Prints reads which are potential chimeras.

These are determined as being soft clipped, not on the edge 
of contigs and at leat 50% mapped with max 2 mismatches.

Input SAM with header.
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

def main():
    contig_len_re = re.compile('SN:(.*)		*LN:([0-9]*)$')
    contig_lens = {}

    # Read the input sam file
    for line in sys.stdin:
        line = line.strip()

        # Check if in header
        if line[0] == '@':
            if line[0:3] == '@SQ':
                contig, contig_len = contig_len_re.findall(line)[0]
                contig_lens[contig] = contig_len
        else:
            # Not in header
            line_split = line.split('\t')
            cigar = line_split[5]
            # Soft clipped at least 20 bases
            if 'S' in cigar:
                soft_clipped, matching = parse_cigar(cigar)
                if soft_clipped > 20:
                    # Matching region should be > 50% of read
                    read_length = len(line_split[9])
                    if matching > read_length/2.0:
                        contig, mapping_start_pos = line_split[2], line_split[7]

                        # mapping should not be on the edges
                        if mapping_start_pos > 300 and int(contig_lens[contig]) - int(mapping_start_pos) > 300:
                            # check nr_mismatches
                            mismatches = int(nr_mismatched_re.findall(line)[0])
                            if int(mismatches) < 2:
                                # check which read 
                                if (int(line_split[1]) / 128) == 1:
                                    read_nr = "2"
                                else:
                                    read_nr = "1"

                                print(line_split[0]+"_" + read_nr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=usage)
    parser.parse_args()
    main()
