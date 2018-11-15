import pandas as pd
import argparse
from collections import defaultdict
import os
import math
from Bio import SeqIO

"""A script to extract rrna sequences for approved MAGs given the metaxa results.

"""

def main(args):
    clustering = pd.read_table(args.clustering_file, sep=',', names=['contig_id', 'cluster_id'], index_col=0)
    taxonomy_df = pd.read_table(args.taxonomy_file, header=None, index_col=0, names=["contig_id", "taxonomy", "identity", "aln_length", "reliability_score"])
    all_approved_prok = pd.read_table(args.all_approved_file, header=None, names=["contig_id"], index_col=0)
    if args.all_approved_euk_file:
        all_approved_euk = pd.read_table(args.all_approved_euk_file, header=None, names=["contig_id"], index_col=0)

    def pair_approved_and_metaxa(all_approved):
        all_approved_set = set(all_approved.index.values)
        contig_to_mag = {}
        taxonomy_df['taxonomy'].fillna('', inplace=True)
        for rrna_contig in taxonomy_df.index.values:
            # Only long contigs are in clustering
            if rrna_contig in clustering.index:
                cluster_id = clustering.loc[rrna_contig]['cluster_id']

                # Only report for approved bins
                if cluster_id in all_approved_set:
                    contig_to_mag[rrna_contig] = cluster_id

        return contig_to_mag
    
    contig_to_mag_prok = pair_approved_and_metaxa(all_approved_prok)

    contig_to_mag_euk = pair_approved_and_metaxa(all_approved_euk)


    def filter_seqs(contig_to_mag, input_fasta ):
        # Keep track if a contig is present several times
        count_per_contig = {}

        output_seqs = []
        # Read prok fasta file
        for seq in SeqIO.parse(input_fasta, "fasta"):
            if seq.id in count_per_contig:
                count_per_contig[seq.id] += 1
            else:
                count_per_contig[seq.id] = 1
            
            # check which MAG it belonged to
            if seq.id in contig_to_mag:
                mag = contig_to_mag[seq.id]
                new_mag_id = args.bin_prefix + str(mag)
                new_seq_id = "{}_{}_{}".format(seq.id, count_per_contig[seq.id], new_mag_id)
                seq.id = new_seq_id

                output_seqs.append(seq)

        return output_seqs

    prok_output_seqs = filter_seqs(contig_to_mag_prok, args.rrna_prok_fasta_file)

    with open(os.path.join(args.outdir, "prok_rRNA.fasta"), 'w') as ofh:
        SeqIO.write(prok_output_seqs, ofh, 'fasta')

    euk_output_seqs = filter_seqs(contig_to_mag_euk, args.rrna_euk_fasta_file)

    with open(os.path.join(args.outdir, "euk_rRNA.fasta"), 'w') as ofh:
        SeqIO.write(euk_output_seqs, ofh, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--clustering_file', help="CONCOCT nocutup merged clustering clustering_nocutup.csv")
    parser.add_argument('--taxonomy_file', help="Metaxa output all_contigs.taxonomy.txt")
    parser.add_argument('--rrna_prok_fasta_file', help="Metaxa output all_contigs.bacteria.fasta and all_contigs.archaea.fasta combined")
    parser.add_argument('--rrna_euk_fasta_file', help="Metaxa output all_contigs.eukaryota.fasta")
    parser.add_argument('--all_approved_file', help="list_of_all_approved_bins_nocutup.tsv")
    parser.add_argument('--all_approved_euk_file', help="list_of_all_approved_bins_nocutup_eukaryotes.tsv")
    parser.add_argument('--bin_prefix', help="Prefix that will be added be fore the integer indicating cluster id")
    parser.add_argument('--outdir', help="A directory for output files")

    args = parser.parse_args()

    main(args)
