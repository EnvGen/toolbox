# coding: utf-8

import pandas as pd
import pysam
import argparse
import sys

def load_mag_read_names(bam_files, mag_names):
    read_sets = {}
    for samfile_path, mag in zip(bam_files, mag_names):
        samfile = pysam.AlignmentFile(samfile_path, 'rb')
        s = set((read.qname, read.is_read1) for read in samfile.fetch())
        read_sets[mag] = s
    return read_sets




def get_aggregated_sets(read_sets, all_mags):
    all_mag_unmapped_reads = None
    atleast_one_mag_unmapped_reads = None
    for mag in all_mags:
        if all_mag_unmapped_reads:
            all_mag_unmapped_reads = all_mag_unmapped_reads.intersection(read_sets[mag])
        else:
            all_mag_unmapped_reads = read_sets[mag]
        if atleast_one_mag_unmapped_reads:
            atleast_one_mag_unmapped_reads = atleast_one_mag_unmapped_reads.union(read_sets[mag])
        else:
            atleast_one_mag_unmapped_reads = read_sets[mag]
    return all_mag_unmapped_reads, atleast_one_mag_unmapped_reads

def main(args):
    mag_read_sets = load_mag_read_names(args.mag_bam_files, args.mag_names)

    output_d = {}
    all_mag_unmapped_reads, atleast_one_mag_unmapped_reads = get_aggregated_sets(mag_read_sets, args.mag_names)
    for samfile_path, sag in zip(args.sag_bam_files, args.sag_names):
        output_d[sag] = {}
        samfile = pysam.AlignmentFile(samfile_path, 'rb')
        s = set((read.qname, read.is_read1) for read in samfile.fetch())

        output_d[sag]["always_missed_and_missed_in_sag"] = len(s.intersection(all_mag_unmapped_reads))

        output_d[sag]["always_missed_not_missed_in_sag"] = len(all_mag_unmapped_reads) - output_d[sag]['always_missed_and_missed_in_sag']

        output_d[sag]["sag_missed_only"] = len(s) - len(s.intersection(atleast_one_mag_unmapped_reads))

        output_d[sag]["atleast_one_not_all_mag_but_missed_in_sag"] =len(s) - output_d[sag]['sag_missed_only'] - output_d[sag]['always_missed_and_missed_in_sag']

        output_d[sag]["atleast_one_not_all_mag_and_not_missed_in_sag"] = len(atleast_one_mag_unmapped_reads) - output_d[sag]['always_missed_and_missed_in_sag'] - output_d[sag]['always_missed_not_missed_in_sag'] - output_d[sag]['atleast_one_not_all_mag_but_missed_in_sag']

    df = pd.DataFrame.from_dict(output_d, orient='index')
    df.to_csv(sys.stdout, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag_bam_files", nargs='*')
    parser.add_argument("--mag_names", nargs='*')
    parser.add_argument("--sag_bam_files", nargs='*')
    parser.add_argument("--sag_names", nargs='*')
    args = parser.parse_args()
    main(args)
