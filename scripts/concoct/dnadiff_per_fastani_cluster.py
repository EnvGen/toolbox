#!/usr/bin/env python
"""
Given a preliminary coarse fastani clustering, run the concoct 
dnadiff script on each cluster respectively.

@author: alneberg
"""
import sys
import os
import argparse
import pandas as pd
from subprocess import call

def main(args):
    # Read fastani clustering
    fastani_df = pd.read_table(args.fastani_clustering, sep=',')
    
    # convert clusters to fasta file paths
    def bin_to_path(bin_id):
        return os.path.join(args.input_dir, "{}.fa".format(bin_id)) 

    fastani_df['filepath'] = fastani_df['bin_id'].apply(bin_to_path)

    # groupby fastani cluster
    for cluster, group_df in fastani_df.groupby('clustering'):
        if len(group_df) == 1:
            print("Cluster {} is a singleton, continuing".format(cluster))
            continue
        # Create dnadiff output folder
        cluster_outdir = os.path.join(args.output_dir, "dnadiff_{}".format(cluster))
        if os.path.isdir(cluster_outdir):
            if args.skip_existing:
                print("Cluster directory {} exists, continuing".format(cluster_outdir))
                continue
        else:
            os.mkdir(cluster_outdir)

        # run dnadiff
        print("Running dnadiff for cluster {}, with {} genomes".format(cluster, len(group_df)))
        fasta_name_file = os.path.join(args.output_dir, 'fasta_names.txt')
        with open(fasta_name_file, 'w') as ofh:
            for fasta_name in list(group_df['bin_id'].values):
                ofh.write(fasta_name + '\n')
        cmd = [args.dnadiff_script, cluster_outdir]
        for filepath in list(group_df['filepath'].values):
            cmd.append(filepath)

        cmd += ["--fasta_names", fasta_name_file, "--skip_plot"]
        if args.skip_dnadiff:
            cmd.append("--skip_dnadiff")
        call(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastani_clustering", help=("Coarse preliminary fastani clustering."))
    parser.add_argument("input_dir", help="Directory where fasta files can be found")
    parser.add_argument("output_dir", help="Directory where dnadiff directories where be placed")
    parser.add_argument("dnadiff_script", help="The script that will be ran")
    parser.add_argument("--skip_dnadiff", action="store_true", help="Skip running dnadiff, assumes it has already been run.")
    parser.add_argument("--skip_existing", action="store_true", help="Skip cluster directories which exists.")
    args = parser.parse_args()

    main(args)
