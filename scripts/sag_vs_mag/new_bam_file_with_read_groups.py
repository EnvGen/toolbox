# coding: utf-8

import pandas as pd
import pysam
import argparse
import sys
import os

def main(args):
    assert not os.path.isfile(args.outfile)

    mag_bam_file = args.mag_bam_file
    mag = args.mag_name
    mag_contig_list = args.mag_contig_list
    sag_bam_file = args.sag_bam_file
    sag = args.sag_name
    sag_unmapped_bam_file = args.sag_unmapped_bam_file


    mag_reads = set()
    long_contig_reads = set()
    metagenome_reads = set()
    date = mag.split('-')[0]

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

    sag_samfile_path = sag_bam_file
    sag_samfile = pysam.AlignmentFile(sag_samfile_path, 'rb')
    sag_reads = set((read.qname, read.is_read1) for read in sag_samfile.fetch())

    sag_unmapped_samfile = pysam.AlignmentFile(sag_unmapped_bam_file, 'rb')
    sag_unmapped_reads = set((read.qname, read.is_read1) for read in sag_unmapped_samfile)


    new_header = sag_samfile.header.copy()

    if not 'RG' in new_header:
        new_header['RG'] = []

    for rg in ['MAG', 'long_metagenome_contig', 'metagenome_contig', 'SAG_only']:
        new_header['RG'].append({'ID': rg})

    outfile = pysam.Samfile(args.outfile, 'wb', header=new_header)
    for read in sag_samfile.fetch():
        read_t = (read.qname, read.is_read1)
        if read_t in mag_reads:
            rg = 'MAG'
        elif read_t in long_contig_reads:
            rg = 'long_metagenome_contig'
        elif read_t in metagenome_reads:
            rg = 'metagenome_contig'
        else:
            rg = 'SAG_only'

        new_tags = []
        for tag in read.tags:
            if tag[0] != 'RG':
                new_tags.append(tag)
        new_tags.append(('RG', rg))
        read.tags = new_tags
        outfile.write(read)
    outfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag_bam_file")
    parser.add_argument("--mag_name")
    parser.add_argument("--mag_contig_list")
    parser.add_argument("--sag_bam_file")
    parser.add_argument("--sag_name")
    parser.add_argument("--sag_unmapped_bam_file")
    parser.add_argument("--outfile")
    args = parser.parse_args()
    main(args)
