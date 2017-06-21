#!/usr/bin/env python
"""A script to sum the values for all genes for each taxonomic annotation."""

import pandas as pd
import argparse
import sys

def output_file_name(output_prefix, level):
    return "{}_{}.tsv".format(output_prefix, level)

def main(args):
    counts_table = pd.read_table(args.counts_table, index_col=0)
    taxonomy = pd.read_table(args.taxonomy_table, index_col=0)

    counts_columns = list(counts_table.columns)
    counts_table['nr_of_genes'] = 1
    counts_columns = ['nr_of_genes'] + counts_columns

    # join the two tables
    joined_table = taxonomy.join(counts_table)

    # For each level, filter away all rows with empty values
    # group by that level and sum all count values.
    level_l = []
    for level in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]:
        level_l.append(level)
        if len(level_l) == 1:
            summed_df = joined_table[joined_table[level].notnull()].groupby(level)[counts_columns].sum()
            summed_df.index.rename("taxonomy")
            summed_df.to_csv(output_file_name(args.output_prefix, level), sep='\t')
        else:
            summed_df = joined_table[joined_table[level].notnull()].groupby(level_l)[counts_columns].sum()
            summed_df["taxonomy"] = [";".join(list(t)) for t in list(summed_df.index)]
            summed_df[["taxonomy"]+counts_columns].to_csv(output_file_name(args.output_prefix, level), sep='\t', index=False)

    joined_table[["superkingdom", "phylum", "class", "order", "family", "genus", "species"]] = joined_table[["superkingdom", "phylum", "class", "order", "family", "genus", "species"]].fillna('')
    joined_table['taxonomy'] = joined_table.superkingdom + ";" + \
                                    joined_table.phylum.map(str) + ";" + \
                                    joined_table['class'].map(str) + ':' + \
                                    joined_table.order.map(str) + ';' + \
                                    joined_table.family.map(str) + ';' + \
                                    joined_table.genus.map(str) +  ';' + \
                                    joined_table.species.map(str)

    summed_df = joined_table.groupby("taxonomy")[counts_columns].sum()
    # For the output file with unique values, add NA to empty cells.
    summed_df.to_csv(output_file_name(args.output_prefix, 'all_levels'), sep='\t', na_rep='NA')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("counts_table")
    parser.add_argument("taxonomy_table")
    parser.add_argument("--output_prefix", default="tax_table", help="Output prefix, default is tax_table")
    args = parser.parse_args()
    main(args)
