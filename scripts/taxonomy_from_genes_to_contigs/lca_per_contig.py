# coding: utf-8
import argparse
import pandas as pd
from skbio import TreeNode
import matplotlib
import numpy as np

from collections import Counter
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def load_ranks(nodes_file):
    ranks = {}
    with open(nodes_file) as nodes_f:
        for line in nodes_f:
            line = line.strip()
            data_fields = line.split('\t')
            ranks[int(data_fields[0])] = {'rank': data_fields[4]}
    return pd.DataFrame.from_dict(ranks, orient='index')


def get_megan_taxassign(megan_result_df, ncbi_tree, ncbi_map, ranks_df, val_cols):
    # Collect information about all ancestors per gene
    megan_taxassign_data = {}
    no_rank = Counter()
    saved_taxa = {}
    default_val_col_values = [None] * len(val_cols)
    for gene_id, row in megan_result_df.iterrows():
        tax_id = row['tax_id']
        data_row = {'tax_id': tax_id}
        data_row.update(zip(val_cols, default_val_col_values))
        if tax_id in saved_taxa:
            data_row.update(saved_taxa[tax_id])
        elif tax_id not in [-2]:
            node = ncbi_tree.find(str(tax_id))

            taxa = ncbi_map['name'].ix[int(node.name)]
            if int(node.name) in ranks_df.index:
                rank = ranks_df['rank'].ix[int(node.name)]
                if rank in val_cols:
                    data_row[rank] = taxa
            else:
                no_rank.update([(taxa, node.name)])

            #Fetch the rank and taxid of all ancestors
            for ancestor in node.ancestors():
                taxa = ncbi_map['name'].ix[int(ancestor.name)]
                if int(ancestor.name) in ranks_df.index:
                    rank = ranks_df['rank'].ix[int(ancestor.name)]
                    if rank in val_cols:
                        data_row[rank] = taxa
                else:
                    no_rank.update([(taxa, ancestor.name)])
            try: saved_taxa[tax_id] = data_row.copy()
            except AttributeError as e: 
                if e.args[0] == "'list' object has no attribute 'copy'": saved_taxa[tax_id] = data_row[:]
            saved_taxa[tax_id].pop('tax_id')
        megan_taxassign_data[gene_id] = data_row
    logging.info("These ranks are missing in the datafile {}".format(no_rank))
    return pd.DataFrame.from_dict(megan_taxassign_data, orient='index')

def get_max(contig_series):
    if len(contig_series):
        val_counts = contig_series.value_counts()
        if len(val_counts):
            return val_counts.index[0], val_counts[0]
        else:
            return None, None
    else:
        return None, None

def partial_contig_lca(df, val_cols):
    try: remaining_cols = val_cols.copy()
    except AttributeError as e: 
        if e.args[0] == "'list' object has no attribute 'copy'": remaining_cols = val_cols[:]
    # Need to check for strange taxa where some levels are skipped
    for col in val_cols[::-1]:
        remaining_cols.remove(col)
        # A mask that finds the special taxa with gaps
        is_null_df = df[df[col].notnull()][remaining_cols].notnull()
        if is_null_df.any().any():
            for gene in is_null_df.index:
                for remaining_col in remaining_cols:
                    if df.ix[gene][remaining_col] is None:
                        df.loc[gene, remaining_col] = "Not Applicable"

    grouped = df.groupby('contig_id')
    nr_assigned_genes_s = grouped[val_cols[0]].agg(lambda x: len(x[x.notnull()]))
    contig_tax_df_temp = pd.DataFrame(index=df['contig_id'].unique())
    contig_tax_df = pd.DataFrame(index=df['contig_id'].unique())
    contig_tax_df_temp['assigned_genes'] = nr_assigned_genes_s
    first = True


    for col in val_cols:
        contig_tax_df_temp[col + '_tuple'] = grouped[col].agg(get_max)
        contig_tax_df_temp[col + '_count'] = contig_tax_df_temp[col + '_tuple'].apply(lambda x: x[1])
        if np.isnan(contig_tax_df_temp[col + '_count'].max()):
            contig_tax_df[col] = None
        else:
            if first:
                to_be_assigned = contig_tax_df_temp[col + '_count'].divide(contig_tax_df_temp['assigned_genes']) >= 0.9
                first = False
            else:
                # to_be_assigned stores if the genes will be assigned on this level or not.
                to_be_assigned = to_be_assigned & (contig_tax_df_temp[col + '_count'].divide(contig_tax_df_temp['assigned_genes']) >= 0.9)
            contig_tax_df[col] = contig_tax_df_temp[col + '_tuple'].apply(lambda x: x[0])[to_be_assigned]

    contig_tax_df[contig_tax_df == 'Not Applicable'] = None
    return contig_tax_df

def apply_contig_assignment(contig_tax_df, df, val_cols):
    classification_lca_df = df['contig_id'].apply(lambda x: contig_tax_df.ix[x])
    classification_lca_df['contig_id'] = df['contig_id']
    return classification_lca_df


def main(args):
    logging.info("Start importing data.")
    df = pd.read_table(args.bed_file, index_col=3, header=None, names=['contig_id', 'start', 'end', 'gene_id'])
    megan_result_df = pd.read_table(args.megan_result, sep=',', index_col=0, header=None, names=['gene_id', 'tax_id'])
    ncbi_tree = TreeNode.read(args.ncbi_tree_file)
    ncbi_map = pd.read_table(args.ncbi_map_file, index_col=0, header=None, names=['id', 'name', 'dunno', 'dunno2'])
    ranks_df = load_ranks(args.nodes_file)

    val_cols = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    logging.info("Finished reading files.")
    megan_taxassign_df = get_megan_taxassign(megan_result_df, ncbi_tree, ncbi_map, ranks_df, val_cols)
    logging.info("Have megan taxassign for each gene")
    # Add megan tax assign to the bed dataframe
    for col in val_cols:
        df[col] = megan_taxassign_df[col]
    no_lca_df = df[ ['contig_id'] + val_cols ]
    no_lca_df.to_csv(args.no_lca_out)

    contig_tax_df = partial_contig_lca(df, val_cols)
    classification_lca_df = apply_contig_assignment(contig_tax_df, df, val_cols)

    logging.info("Finished partial lca")

    lca_output_df = classification_lca_df[['contig_id'] + val_cols]
    lca_output_df.to_csv(args.lca_out)


def partial_lca(df, val_cols):
    # Partial LCA Classification for all genes within a contig
    return_cols = ['contig_id'] + val_cols
    try: classification_lca_df = df[return_cols].copy()
    except AttributeError as e: 
        if e.args[0] == "'list' object has no attribute 'copy'": classification_lca_df = df[return_cols][:]
    classification_lca_df.loc[df.index][return_cols] = np.NaN
    for col in val_cols:
        df.groupby(['contig_id', col]).size()
    for contig_id, contig_df in df.groupby('contig_id'):
        assigned_genes = contig_df[np.any(contig_df[val_cols].notnull(), axis=1)]
        time_to_abort = False
        for i, col in enumerate(val_cols):
            if time_to_abort:
                break

            assigned_gene_counts = assigned_genes.groupby(col).size()
            if len(assigned_gene_counts) > 0:
                argmax = assigned_gene_counts.argmax()
                max_count = assigned_gene_counts[argmax]

                if max_count/float(len(assigned_genes)) > 0.9:
                    assigned_taxa = argmax
                    classification_lca_df.loc[contig_df.index][col] = assigned_taxa
                else:
                    time_to_abort = True
            else:
                time_to_abort = True
    return classification_lca_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--ncbi_tree_file", help="ncbi.tre")
    parser.add_argument("--ncbi_map_file", help="ncbi.map")
    parser.add_argument("--nodes_file", help="nodes.dmp")
    parser.add_argument("--megan_result", help="Megan output as csv with gene_id and tax_id columns")
    parser.add_argument("--bed_file")
    parser.add_argument("--no_lca_out")
    parser.add_argument("--lca_out")
    args = parser.parse_args()

    main(args)
