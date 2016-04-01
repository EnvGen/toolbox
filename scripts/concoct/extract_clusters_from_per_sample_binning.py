#!/usr/bin/env python
"""extract_clusters_from_per_sample_binning.py 

create a directory structure where each directory contains links to 
all fasta files within a certain cluster. If --special_fastas are given, 
these are used to name the output directory.
"""

import argparse
import sys
import os
import pandas as pd

def main(args):
    df = pd.read_csv(args.cluster_file, header=None, names=["fasta_file", "cluster_id"], sep='\t')
    
    # Annotate the special fastas

    # Check that no cluster have more than one type of special fastas

    # Check that no special fasta type is split over more than one cluster

    cluster_to_contigs = defaultdict(list)
    for i, row in df.iterrows():
        cluster_to_contigs[row['cluster_id']].append(row['contig_id'])
    
    for cluster_id, contig_ids in cluster_to_contigs.iteritems():
        output_file = os.path.join(args.output_path, "{0}.fa".format(cluster_id))
        seqs = [all_seqs[contig_id] for contig_id in contig_ids] 
        with open(output_file, 'w') as ofh:
            SeqIO.write(seqs, ofh, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", help="Top output directory where cluster directories will be created.")
    parser.add_argument("cluster_file", help="Concoct output cluster file")
    parser.add_argument("--special_files", help=("A file containing the names of the fasta " 
            "files within the clustering that is considered 'special'."
            " These will be used to name the output directories for each cluster."))
    args = parser.parse_args()

    main(args)
