import pandas as pd
import argparse
from collections import defaultdict
import os
import math

def main(args):
    clustering = pd.read_table(args.clustering_file, sep=',', names=['contig_id', 'cluster_id'], index_col=0)
    taxonomy_df = pd.read_table(args.taxonomy_file, header=None, index_col=0, names=["contig_id", "taxonomy", "identity", "aln_length", "reliability_score"])
    all_approved_prok = pd.read_table(args.all_approved_file, header=None, names=["contig_id"], index_col=0)
    if args.all_approved_euk_file:
        all_approved_euk = pd.read_table(args.all_approved_euk_file, header=None, names=["contig_id"], index_col=0)

    def pair_approved_and_metaxa(all_approved):
        all_approved_set = set(all_approved.index.values)
        unapproved_rrna = defaultdict(int)
        all_clusters_found = set()
        approved_rrna = []
        taxonomy_df['taxonomy'].fillna('', inplace=True)
        for rrna_contig in taxonomy_df.index.values:
            if rrna_contig in clustering.index:
                cluster_id = clustering.loc[rrna_contig]['cluster_id']
                metaxa_val = taxonomy_df.loc[rrna_contig]
                if cluster_id in all_approved_set:
                    tax_dict = metaxa_val.to_dict()
                    tax_dict['cluster_id'] = cluster_id
                    tax_dict['contig_id'] = rrna_contig
                    approved_rrna.append(tax_dict)
                    all_clusters_found.add(cluster_id)

        for cluster_id in all_approved_set:
            if cluster_id not in all_clusters_found:
                approved_rrna.append({'cluster_id': cluster_id, 'contig_id': '', 'taxonomy': '', 'identity': '', 'aln_length': '', 'reliability_score': ''})

        return approved_rrna
    
    approved_rrna_prok = pair_approved_and_metaxa(all_approved_prok)
    approved_stats_prok_df = pd.DataFrame(approved_rrna_prok)

    approved_rrna_euk = pair_approved_and_metaxa(all_approved_euk)
    approved_stats_euk_df = pd.DataFrame(approved_rrna_euk)

    columns = ["cluster_id", "contig_id", "taxonomy", "identity", "aln_length", "reliability_score"]
    with open(os.path.join(args.outdir, 'approved_prok.tsv'), 'w') as ofh:
        if len(approved_stats_prok_df):
            approved_stats_prok_df[columns].to_csv(ofh, sep='\t', index=False)
        else:
            ofh.write('\t'.join(columns) + '\n')
    with open(os.path.join(args.outdir, 'approved_euk.tsv'), 'w') as ofh:
        if len(approved_stats_euk_df):
            approved_stats_euk_df[columns].to_csv(ofh, sep='\t', index=False)
        else:
            ofh.write('\t'.join(columns) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--clustering_file', help="CONCOCT nocutup merged clustering clustering_nocutup.csv")
    parser.add_argument('--taxonomy_file', help="Metaxa output all_contigs.taxonomy.txt")
    parser.add_argument('--all_approved_file', help="list_of_all_approved_bins_nocutup.tsv")
    parser.add_argument('--all_approved_euk_file', help="list_of_all_approved_bins_nocutup_eukaryotes.tsv")
    parser.add_argument('--outdir', help="A directory for output files")

    args = parser.parse_args()

    main(args)
