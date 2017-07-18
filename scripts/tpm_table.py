#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
"""A script to calculate TPM values for contigs or genes based on count files

TPM values are defined as in Wagner et al (Theory in Biosciences) 2012. 

      rg x rl x 10^6
TPM = --------------
        flg x T

rg: reads mapped to gene g
rl: read length
flg: feature length 
T: sum of rgxrl/flg for all genes
"""
import sys, pandas as pd, argparse, logging
import re
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def gene_lengths_from_gff(gff_file):
    gene_id_regex = re.compile('ID=([a-zA-Z_\-0-9]*);')
    gene_lengths = {}
    with open(gff_file) as fh:
        for line in fh:
            gene_id = gene_id_regex.findall(line)[0]
            gene_lengths[gene_id] = abs(int(line.split('	')[4]) - int(line.split('	')[3])) + 1
    return pd.Series(gene_lengths)

def main(args):
    logging.info("Reading sample info")
    sample_info = pd.read_table(args.sample_info, header=None, index_col=0, names=['avg_read_len'])
    logging.info("Reading gene lengths from gff")

    gene_lengths = gene_lengths_from_gff(args.gff)

    df = pd.DataFrame()
    first = True
    for fn, sample_name in zip(args.coverage_files, args.sample_names):
        logging.info("Calculating TPM for "+ sample_name)
        ## Read counts per gene for sample
        rg = pd.read_table(fn, index_col=0, header=None, names=['gene_id', 'count'], compression=args.input_compression)
        ## Intersect with genes in the gene length file
        rg = rg.loc[list(set(gene_lengths.index).intersection(set(rg.index)))]
        gene_lengths = gene_lengths.loc[list(rg.index)]
        ## Average read length for sample
        rl = sample_info.ix[sample_name,'avg_read_len']
        ## Calculate T for sample
        T = rl * rg['count'].divide(gene_lengths).sum()
        ## Calculate TPM for sample
        tpm = ((1e6*rl)/float(T))*(rg['count'].divide(gene_lengths))
        ## Create dataframe
        TPM = pd.DataFrame(tpm,columns=[sample_name])

        ## Add gene length as the first column
        if first:
            first = False
            df['gene_length'] = gene_lengths
        
        ## Concatenate to results
        df = pd.concat([df,TPM],axis=1)
    df.index.name = 'gene_id'
    ## Write to file
    df.to_csv(sys.stdout, sep='\t')
    logging.info("Done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--sample_names', nargs='*', 
            help="Sample names, in the same order as coverage_files")
    parser.add_argument('-c', '--coverage_files', nargs='*', 
            help="Coverage files with tab separated values: 'sequence id, count'")
    parser.add_argument('--gff',
            help=("GFF version 2 file"))
    parser.add_argument('-i', '--sample_info', 
            help="Tab separated values 'sample_id', 'avg_read_length'")
    parser.add_argument('-l', '--gene_lengths',
            help="Gene lengths in a tsv file")
    parser.add_argument("--input_compression", default=None, choices=[None, 'gzip'], 
            help="Compression type for input coverage files. Default=None, use 'gzip', for gzipped files.")
    args = parser.parse_args()
    main(args)
