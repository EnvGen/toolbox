# coding: utf-8
import pandas as pd
from collections import defaultdict

def main(args):
    clustering = pd.read_table(args.clustering_file, sep=',', names=['contig_id', 'cluster_id'], index_col=0)
    taxonomy_df = pd.read_table(args.taxonomy_file, header=None, index_col=0, names=["contig_id", "taxonomy", "bla", "bla1", "bla2"])
    all_approved = pd.read_table(args.all_approved_file, header=None, names=["contig_id"], index_col=0)
    checkm_taxonomy = pd.read_table(args.checkm_taxonomy_file, index_col=0)

    all_approved_set = set(all_approved.index.values)
    unapproved_rrna = defaultdict(int)
    approved_rrna = {}
    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for rrna_contig in taxonomy_df.index.values:
        if rrna_contig in clustering.index:
            cluster_id = clustering.loc[rrna_contig]['cluster_id']
            if cluster_id in all_approved_set:
                checkm_val = checkm_taxonomy.loc[cluster_id]['Taxonomy'].split(';')
                metaxa_val = taxonomy_df.loc[rrna_contig]['taxonomy'].split(';')

                metaxa_val = fix_strange_metaxa_vals(metaxa_val)


                matched_level = None
                for i, level in enumerate(levels):
                    checkm_level_val, metaxa_level_val = None, None
                    if len(checkm_val) > i and len(metaxa_val) > i:
                        checkm_level_val = checkm_val[i][3:]
                        metaxa_level_val = metaxa_val[i]

                        if level == 'species':
                            metaxa_level_val = metaxa_val[i].replace(' ', '_')
                        if checkm_level_val == metaxa_level_val:
                            matched_level = i
                        else:
                            break
                    else:
                        matched_level = i-1
                        break
                if cluster_id not in approved_rrna:
                    approved_rrna[cluster_id] = {'matching': 0, 'not matching': 0}

                if matched_level >= 3:
                    approved_rrna[cluster_id]['matching'] += 1
                else:
                    approved_rrna[cluster_id]['not matching'] += 1
                #print(most_detailed_level_checkm, most_detailed_level_metaxa)
                #print(most_detailed_matched_level)
                #print(taxonomy_df.loc[rrna_contig]['taxonomy'], checkm_taxonomy.loc[cluster_id]['Taxonomy'])
            else:
                unapproved_rrna[cluster_id] += 1

    for cluster_id in all_approved_set:
        if cluster_id not in approved_rrna:
            approved_rrna[cluster_id] = {'matching': 0, 'not matching': 0}

    approved_stats_df = pd.DataFrame.from_dict(approved_rrna, orient='index')

    unapproved_stats_df = pd.DataFrame.from_dict(unapproved_rrna, orient='index')
    unapproved_stats_df.columns = ['nr_rrna']

    print(approved_stats_df)
    print(unapproved_stats_df)

    # Things to output:
    #
    #     Number of approved genomes with matching rrna
    #     Number of approved genomes with unmatching rrna
    #     Number of rrna genes per bin
    #
    #     Number of approved genomes with > 0 matching rrna and < 2 unmatching rrna
    #     Matching is counted at order level
    #


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--approved_bin_dir', help="Directory where bins can be found, will be copied to directories decided on taxonomy.")
    parser.add_argument('--all_approved_file', help="e.g. ../../Data/test_binning_and_16s_combo/list_of_all_approved_bins_nocutup.tsv")
    parser.add_argument('--checkm_taxonomy_file', help="e.g. ../../Data/test_binning_and_16s_combo/checkm_tree_qa.tsv")

    args = parser.parse_args()

    main(args)
