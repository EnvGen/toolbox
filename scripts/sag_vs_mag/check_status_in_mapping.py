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

left_clipped_re = re.compile('^([0-9]+)S')

nr_mismatched_re = re.compile('XM:i:([0-9]+)')

def parse_cigar(cigar):
    soft_clipped_matches = soft_clip_re.findall(cigar)
    matched_matches = matched_re.findall(cigar)    
    soft_clipped_bases = sum([int(x) for x in soft_clipped_matches])
    matched_bases = sum([int(x) for x in matched_matches])

    left_clipped_matches = left_clipped_re.findall(cigar)
    left_clipped_bases = sum([int(x) for x in left_clipped_matches])

    return soft_clipped_bases, matched_bases, left_clipped_bases


"""
chimeric or not
soft clipped (> 20S) or not 
overlapping edge or not
mapped or not

              soft clipped   not soft clipped    overlapping edge    not overlapping edge   overlapping edge and soft clipped    mapped    not mapped
    chimeric: 
Non-chimeric:  

"""

def print_stats(name, soft_clipped_set, overlapping_set, mapped_set, complete_set):
    "soft clipped   not soft clipped    overlapping edge    not overlapping edge   overlapping edge and soft clipped    mapped    not mapped"

    print("{},{},{},{},{},{},{},{}".format(name, len(soft_clipped_set), len(complete_set - soft_clipped_set), len(overlapping_set), len(complete_set - overlapping_set), len(overlapping_set & soft_clipped_set), len(mapped_set), len(complete_set - mapped_set)))


def main(args):
    chim_reads = set()
    for line in open(args.subset_reads, 'r'):
        read = line.strip()
        chim_reads.add(read)

    mapped_reads = set()
    soft_clipped_reads = set()
    overlapping_edge_reads = set()
    all_reads = set()
    contig_lens = {}

    contig_len_re = re.compile('SN:(.*)		*LN:([0-9]*)$')

    # Read the input sam files
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
            if (int(line_split[1]) / 128) == 1:
                read_id = line_split[0] + "_2"
            else:
                read_id = line_split[0] + "_1"

            all_reads.add(read_id)
            
            mapping_bit = (((((int(line_split[1]) % 128) % 64) % 32) % 16) % 8) / 4
            if mapping_bit != 1:
                mapped_reads.add(read_id)
                cigar = line_split[5]

                soft_clipped, matching, left_clipped_bases = parse_cigar(cigar)
                if soft_clipped > 20:
                    soft_clipped_reads.add(read_id)

                read_length = len(line_split[9])
                 
                contig, mapping_start_pos = line_split[2], line_split[7]
                
                if int(mapping_start_pos) - left_clipped_bases < 0:
                    overlapping_edge_reads.add(read_id)
		   
                if int(mapping_start_pos) + read_length > contig_lens[contig]:
                    overlapping_edge_reads.add(read_id)

    # Aim to print> 
    # "soft clipped   not soft clipped    overlapping edge    not overlapping edge   overlapping edge and soft clipped    mapped    not mapped"
    print("name,soft clipped, not soft clipped,overlapping edge,not overlapping edge,overlapping edge and soft clipped,mapped,not mapped")
    print_stats("chimeric_reads", soft_clipped_reads & chim_reads, overlapping_edge_reads & chim_reads, mapped_reads & chim_reads, chim_reads)
    print_stats("other_reads", soft_clipped_reads - chim_reads, overlapping_edge_reads - chim_reads, mapped_reads - chim_reads, all_reads - chim_reads) 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('subset_reads')
    args = parser.parse_args()

    main(args)
