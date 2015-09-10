#!/usr/bin/env python

"""Combine_classifications.py: fills in taxonomic classification from separate tsv.

From tab-separated blast files corresponding to fwd and rev reads,
the 5% top HSP making cutoffs for evalue, coverage and identity 
in both files are used for LCA determination"""

from collections import defaultdict
import argparse
import csv
import re


__author__ = "Yue O Hu and Luisa W Hugerth"
__email__ = "luisa.hugerth@scilifelab.se"

def parsetax(taxdict, leveldict, taxfile, level):
	with open(taxfile) as csvfile:
		reader = csv.reader(csvfile, delimiter="\t")
		for row in reader:
			query = row[0]
			tax = row[1]		
			if query not in taxdict or taxdict[query] == "Unclassified":
				taxdict[query] = tax
				leveldict[query] = level
	return taxdict, leveldict


def main(infiles, names):
	filelist = infiles.split(",")
	namelist = names.split(",")
	count = 0
	taxdict = dict()
	leveldict = dict()
	for infile in filelist:
		taxdict, leveldict = parsetax(taxdict, leveldict, infile, namelist[count])
		count += 1
	for query, tax in taxdict.iteritems():
		print query + "\t" + leveldict[query] + "\t" + tax 

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Combines taxonomic assignment at different stringency levels')
	parser.add_argument('-i', '--infiles', help='Paths to relevant tsv taxonomies (outputs of taxonomy_blast_parser.py) in priority order,\
													\ separated by ","')
	parser.add_argument('-n', '--names', help='Short names for the taxonomines, in the same order, separated by ","')
	args = parser.parse_args()

	main(args.infiles, args.names)

