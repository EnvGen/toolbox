#!/usr/bin/env python
"""construct_ena_sequencing_runs_table.py 

Based on a folder with uploaded files and a template, construct table ready to be submitted.
This script is not general but is very niched to the NGI/Uppmax scenario.
Ideally the user is expected to copy this script and edit it to suit the users needs.
"""

import argparse
import sys
import os
import glob
from os.path import join as opj
import pandas as pd
from collections import defaultdict
import gzip

# Need to fetch file name and link it to sample id
# Fetch the md5sum already calculated

def main(args):
    md5sum_df = pd.read_table(args.md5_summary, sep=',', header=None, names=['file_name', 'md5sum'], index_col=0)

    insert_sizes = pd.read_table(args.insert_size, index_col=0)
    info_d = {}

    for sample_dir in args.sample_dirs:
        for R1_run_file in glob.glob(opj(sample_dir, "*", "*_R1*.fastq.gz")):
            R1_file_name=os.path.basename(R1_run_file)
            sample_name="_".join(R1_file_name.split('_')[0:2])
            run_id="_".join(R1_file_name.split('_')[0:4])
            run_info = {}
            run_info['sample_accession'] = sample_name
            run_info['library_name'] = run_id
            is_series = insert_sizes.ix[sample_name]['Avg. FS']
            try:
                run_info['insert_size'] = int(round(is_series))
            except TypeError:
                run_info['insert_size'] = int(round(insert_sizes[insert_sizes['Lib QC'] == 'PASSED'].ix[sample_name]['Avg. FS']))

            run_info['forward_file_name'] = R1_file_name
            run_info['forward_file_md5'] = md5sum_df.loc[R1_file_name]['md5sum']
            R2_file_name = R1_file_name.replace("R1", "R2")
            run_info['reverse_file_name'] = R2_file_name
            run_info['reverse_file_md5'] = md5sum_df.loc[R2_file_name]['md5sum']
            run_info['library_source'] = 'METAGENOMIC'
            run_info['library_selection'] = 'RANDOM'
            run_info['library_strategy'] = 'WGS'
            run_info['library_construction_protocol'] = 'Rubicon Thruplex'
            run_info['instrument_model'] = 'Illumina HiSeq 2500'
            run_info['file_type'] = 'fastq'
            run_info['library_layout'] = 'PAIRED'

            
            info_d[run_id] = run_info 
    all_columns_sorted = ['sample_accession', 'library_name', 'library_source', 'insert_size', \
            'library_selection', 'library_strategy', 'library_construction_protocol', 'instrument_model', \
            'file_type', 'library_layout', 'insert_size', \
            'forward_file_name', 'forward_file_md5', 'reverse_file_name', 'reverse_file_md5']

    df = pd.DataFrame.from_dict(info_d, orient='index')
    df[all_columns_sorted].to_csv(sys.stdout, index=False, sep='\t', header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("md5_summary", help="Table (csv) with md5sum values for all files")
    parser.add_argument("insert_size", help="Table with insert size per sample")
    parser.add_argument("sample_dirs", nargs='*', help="Directories where read files are located in subdirs")
    args = parser.parse_args()

    main(args)
