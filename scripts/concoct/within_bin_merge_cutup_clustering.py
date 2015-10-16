#!/usr/bin/env python
"""
With contigs cutup with cut_up_fasta.py as input, sees to that the consequtive
parts of the original contigs are merged.

prints result to stdout.

@author: alneberg
"""
from __future__ import print_function
import sys
import os
import argparse
from Bio import SeqIO

def original_contig_name(s):
    """Transform s to the original contig name"""
    n = s.split(".")[-1]
    try:
        int(n)
    except:
        return s, 0
    # Only small integers are likely to be 
    # indicating a cutup part.
    if int(n) < 1000:
        return ".".join(s.split(".")[:-1]), n
    else:
        # A large n indicates that the integer
        # was part of the original contig
        return s, 0


def original_contig_name_special(s):
    """Transform s to the original contig name according to the special Baltic Mag paper"""
    n = s.split(".")[-1].split('_')[0]
    try:
        int(n)
    except:
        return s, 0, None
    # Only small integers are likely to be 
    # indicating a cutup part.
    if int(n) < 1000:

        return ".".join(s.split(".")[:-1]), int(n), "_".join(s.split(".")[-1].split('_')[1:])
    else:
        # A large n indicates that the integer
        # was part of the original contig
        return s, 0, None

def main(args):

    all_seqs = {}
    for i, seq in enumerate(SeqIO.parse(args.fasta_file, "fasta")):
        all_seqs[seq.id] = seq

    first_iteration = True
    merged_contigs_stack = []

    for seq_id in sorted(all_seqs.keys()):
        original_contig_id, part_id, left_over = original_contig_name_special(seq_id)
 
        if first_iteration:
            previous_orig_id = original_contig_id
            previous_part_id = part_id
            merged_contig = [all_seqs[seq_id]]
            first_iteration = False
            continue
        if previous_orig_id == original_contig_id:
            # This assumes they are sorted
            if part_id - 1 == previous_part_id:
                # Consequtive parts can be merged
                # Add to previous group
                merged_contig.append(all_seqs[seq_id])
            else:
                # Add finished group to stack
                merged_contigs_stack.append(merged_contig)

                # Create a new group
                merged_contig = [all_seqs[seq_id]]
        else:
            merged_contigs_stack.append(merged_contig)
            merged_contig = [all_seqs[seq_id]]

        previous_orig_id = original_contig_id 
        previous_part_id = part_id

    if merged_contig:
        merged_contigs_stack.append(merged_contig)

    for merged_contig in merged_contigs_stack:
        if len(merged_contig) > 1:
            first_contig = merged_contig[0]
            last_contig = merged_contig[-1]

            original_contig_id, first_part_id, left_over = original_contig_name_special(first_contig.id)
            original_contig_id, last_part_id, left_over  = original_contig_name_special(last_contig.id)

            new_id = original_contig_id + ".{0}_to_{1}_{2}".format(first_part_id, last_part_id, left_over)

            first_contig.id = new_id
            first_contig.description = new_id

            for contig_part in merged_contig[1:]:
                first_contig.seq += contig_part.seq
            SeqIO.write([first_contig], sys.stdout, 'fasta')
        else:
            SeqIO.write(merged_contig, sys.stdout, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_file", help=("Input Fasta file."))
    args = parser.parse_args()

    main(args)
