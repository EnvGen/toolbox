#!/usr/bin/env python
usage="""Checks the status for a subset of reads within a bam file.

Three categories:
    Not Mapped
    Mapping to < 1000 bp contig
    Mapping to > 1000 bp contig (not MAG)
    Mapping to MAG contig
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

    MAG_contigs = set()
    for line in open(args.mag_contigs, 'r'):
        contig = line.strip()
        MAG_contigs.add(contig)

    not_mapped = set()
    short_mapped = set()
    long_mapped = set()
    MAG_mapped = set()

    contig_len_re = re.compile('SN:(.*)		*LN:([0-9]*)$')
    contig_lens = {}

    # Read the input sam file
    for line in sys.stdin:
        line = line.strip()

        # Check if in header
        if line[0] == '@':
            if line[0:3] == '@SQ':
                contig, contig_len = contig_len_re.findall(line)[0]
                contig = contig.split('_')[0]
                contig_lens[contig] = int(contig_len)
        else:
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
                    contig = line_split[2]
                    contig = contig.split('_')[0]
                    if contig in MAG_contigs:
                        MAG_mapped.add(read_id)
                    elif contig_lens[contig] >= 1000:
                        long_mapped.add(read_id)
                    else:
                        short_mapped.add(read_id)

    print("Not Mapped\tShort Mapped\tLong Mapped\tMAG mapped")
    print("{}\t{}\t{}\t{}".format(len(not_mapped), len(short_mapped),
        len(long_mapped), len(MAG_mapped)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('subset_reads')
    parser.add_argument('mag_contigs')
    args = parser.parse_args()

    main(args)
