#!/usr/bin/env python
usage="""Filter potential chimeric reads to remove correctly none-clipped. 
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
    for line in open(args.potential_chimera, 'r'):
        read = line.strip()
        all_reads.add(read)

    # Read the input sam files
    for line in sys.stdin:
        line = line.strip()

        # Check if in header
        if line[0] != '@':
            # Not in header
            
            # Ignore reads which are note potential chimeras
            # Read is most likely not chimeric
            line_split = line.split('\t')
            if (int(line_split[1]) / 128) == 1:
                read_id = line_split[0] + "_2"
            else:
                read_id = line_split[0] + "_1"

            if read_id in all_reads:
                cigar = line_split[5]
                # Matching min 95% of bases
                soft_clipped, matching = parse_cigar(cigar)
                    
                read_length = len(line_split[9])
                mismatches = int(nr_mismatched_re.findall(line)[0])
                
                if (matching > read_length*0.95) and (int(mismatches) < 5):
                    all_reads.remove(read_id)

    print("\n".join(list(all_reads)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('potential_chimera')
    args = parser.parse_args()

    main(args)
