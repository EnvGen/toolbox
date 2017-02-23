#!/usr/bin/env python
"""A script to quantify duplicated segments

"""

import sys
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

def read_gff_for_zero_cov_contigs(gff_file, df):
    data = {}
    covered_contigs = set(df['contig_id'].unique())
    with open(gff_file, 'r') as fh:
        for i, line in enumerate(fh):
            if i!=0 and line.startswith('#'):
                line = line.strip()
                _, contig, start, end = line.split(' ')
                if not (contig in covered_contigs):
                    data[contig] = (start,end)
            elif line.startswith('#'):
                continue
            else:
                break
    return data

def values_from_original_bed(bed_df, values, sample_name, first=True):
    coverage_name = "coverage_{}".format(sample_name)
    for ix, row in bed_df.iterrows():
        contig_len = int(row['contig_id'].split('_')[3])
        segment_len = row['end'] - row['start']
        keys = zip([row["contig_id"]]*segment_len, list(range(row['start'], row['end'])))
        if first:
            vals = [{coverage_name: row['coverage'], "duplicated": 0, "distance_from_end": contig_len - x -1, "distance_from_start": x} for x in range(row['start'], row['end'])]
            values.update(dict(zip(keys,vals)))
        else:
            for key in keys:
                values[key][coverage_name] = row['coverage']
    return values

def values_from_original_bed_speedup(bed_df, values, sample_name, first=True):
    coverage_name = "coverage_{}".format(sample_name)
    coverage = []
    duplicated = []
    distance_from_end = []
    distance_from_start = []
    for ix, row in bed_df.iterrows():
        contig_len = int(row['contig_id'].split('_')[3])
        segment_len = row['end'] - row['start']
        keys = zip([row["contig_id"]]*segment_len, list(range(row['start'], row['end'])))
        if first:
            vals = [{coverage_name: row['coverage'], "duplicated": 0, "distance_from_end": contig_len - x -1, "distance_from_start": x} for x in range(row['start'], row['end'])]
            values.update(dict(zip(keys,vals)))
        else:
            for key in keys:
                values[key][coverage_name] = row['coverage']
    return values


def values_from_zero_cov_contigs(zero_cov_contigs, values_all, sample_name, first=True):
    coverage_name = "coverage_{}".format(sample_name)
    for contig, start_end in zero_cov_contigs.iteritems():
        contig_len = int(start_end[1]) - int(start_end[0]) + 1
        keys = zip([contig]*contig_len, list(range(contig_len)))
        if first:
            vals = [{coverage_name: 0, "duplicated": 0, "distance_from_end": contig_len - x -1, "distance_from_start": x} for x in range(contig_len)]
            values_all.update(dict(zip(keys,vals)))
        else:
            for key in keys:
                values_all[key][coverage_name] = 0
    return values_all

def update_values_with_duplication(dup_regions_df, values_all):
    for ix, row in dup_regions_df.iterrows():
        contig_len = int(row['contig_id'].split('_')[3])
        segment_len = row['end'] - row['start']

        keys = zip([row['contig_id']]*segment_len, list(range(segment_len)))
        for key in keys:
            values_all[key]['duplicated'] += 1
    return values_all

def main(args):
    values = {}
    sag = args.sag
    for i, mag_date in enumerate(args.mag_dates):
        df = pd.read_table(args.coverage_files[i], names=['contig_id', 'start', 'end','coverage'])
        values = values_from_original_bed(df, values, mag_date, not i)
        if not i:
            zero_cov_contigs = read_gff_for_zero_cov_contigs(args.gff_file, df)
            values = values_from_zero_cov_contigs(zero_cov_contigs, values, mag_date, not i)

    dup_regions_df = pd.read_table(args.sag_duplicates)
    values = update_values_with_duplication(dup_regions_df, values)

    base_cov_df = pd.DataFrame.from_dict(values, orient='index')
    base_cov_df.fillna(0, inplace=True)
    base_cov_df.index.names = ['contig_id', 'position']
    duplicates = base_cov_df['duplicated'] > 0

    middle_bases = (base_cov_df['distance_from_end'] > 200) & (base_cov_df['distance_from_start'] > 200)

    outliers_d = defaultdict(dict)
    for mag_date in args.mag_dates:
        coverage_name = "coverage_{}".format(mag_date)
        outliers_d[coverage_name][True] = base_cov_df[coverage_name][duplicates] > base_cov_df[coverage_name][duplicates].median()+base_cov_df[coverage_name][duplicates].std()*3
        outliers_d[coverage_name][False] = base_cov_df[coverage_name][~duplicates] > base_cov_df[coverage_name][~duplicates].median()+base_cov_df[coverage_name][~duplicates].std()*3

    results = {}

    for mag_date in args.mag_dates:
        coverage_name = "coverage_{}".format(mag_date)
        results[(sag, mag_date)] = {}
        outliers = outliers_d[coverage_name][False]
        results[(sag, mag_date)]["Not duplicated, not outliers"] = base_cov_df[middle_bases & (~duplicates) & (~outliers)][coverage_name].median()

        outliers = outliers_d[coverage_name][False]
        results[(sag, mag_date)]["Not duplicated, all values"] = base_cov_df[middle_bases & (~duplicates)][coverage_name].median()

        outliers = outliers_d[coverage_name][True]
        results[(sag, mag_date)]["Duplicated, not outliers"] = base_cov_df[middle_bases & (duplicates) & (~outliers)][coverage_name].median()
        results[(sag, mag_date)]["Duplicated, nr bases, not outliers"] = len(base_cov_df[middle_bases & (duplicates) & (~outliers)])

        outliers = outliers_d[coverage_name][True]
        results[(sag, mag_date)]["Duplicated, all values"] = base_cov_df[middle_bases & (duplicates)][coverage_name].median()

        results[(sag, mag_date)]["Duplicated, nr bases, all values"] = len(base_cov_df[middle_bases & (duplicates)])


    df = pd.DataFrame.from_dict(results, orient='index')
    df.to_csv(sys.stdout, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sag")
    parser.add_argument("--mag_dates", nargs='*')
    parser.add_argument("--coverage_files", nargs='*')
    parser.add_argument("--sag_duplicates")
    parser.add_argument("--gff_file")

    args = parser.parse_args()
    main(args)

