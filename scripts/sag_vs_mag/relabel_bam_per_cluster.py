# coding: utf-8

import numpy as np
import pandas as pd
import pysam
import argparse
import sys

# Group SAG reads into groups:
# 1. Reads that map to both SAG and at least one MAG
# 2. Reads that map to at least one metagenome contig shorter than 1000 but never a MAG contig and map to SAG
# 3. Reads that map to metagenome but never a contig longer than 1000, and map to SAG
# 4. Reads that do not map to metagenome but to SAG x

def main(args):
    arg_lists = [args.mag_bam_files, args.mag_names, args.mag_contig_lists]

    for mag_bam_file, mag, mag_contig_list in zip(*arg_lists):
        mag_reads = set()
        long_contig_reads = set()
        metagenome_reads = set()

        samfile = pysam.AlignmentFile(mag_bam_file, 'rb')
        bin_contigs = pd.read_table(mag_contig_list, header=None, names=["contig_id"])
        all_mag_contigs = set(bin_contigs['contig_id'].values)
        for contig, contig_len in zip(samfile.references, samfile.lengths):
            contig_read_set = set((read.qname, read.is_read1) for read in samfile.fetch(contig))
            metagenome_reads |= contig_read_set
            if contig_len >= 1000:
                long_contig_reads |= contig_read_set
            if contig in all_mag_contigs:
                mag_reads |= contig_read_set
    
    long_contig_reads -= mag_reads
    metagenome_reads -= mag_reads
    metagenome_reads -= long_contig_reads

    sag_samfile = pysam.AlignmentFile(args.sag_bam_file, 'rb')
    header = sag_samfile.header
    header['RG'].append({"ID": "MAG_mapped"})
    header['RG'].append({"ID": "Long_non_MAG_contig_mapped"})
    header['RG'].append({"ID": "Short_metagenome_contig_mapped"})
    header['RG'].append({"ID": "Metagenome_unmapped"})
    with pysam.AlignmentFile(args.sag_new_bam, "wb", header=header) as outf:
        for read in sag_samfile.fetch():
            old_tags = read.tags
            if (read.qname, read.is_read1) in metagenome_reads:
                new_tags = old_tags + [("RG", "Short_metagenome_contig_mapped")]
            elif (read.qname, read.is_read1) in mag_reads:
                new_tags = old_tags + [("RG", "MAG_mapped")]
            elif (read.qname, read.is_read1) in long_contig_reads:
                new_tags = old_tags + [("RG", "Long_non_MAG_contig_mapped")]
            else:
                new_tags = old_tags + [("RG", "Metagenome_unmapped")]
            read.tags = new_tags
            outf.write(read)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag_bam_files", nargs='*')
    parser.add_argument("--mag_names", nargs='*')
    parser.add_argument("--mag_contig_lists", nargs='*')
    parser.add_argument("--sag_bam_file")
    parser.add_argument("--sag_new_bam")
    args = parser.parse_args()
    main(args)
