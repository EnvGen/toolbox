#!/usr/bin/env python

"""Taxonomy_blast_parser.py: finds the best taxonomic classification from 2 blast files.

From tab-separated blast files corresponding to fwd and rev reads,
the 5% top HSP making cutoffs for evalue, coverage and identity 
in both files are used for LCA determination"""

from __future__ import division
from collections import defaultdict
from os.path import commonprefix
import argparse
import csv
import math
import re


__author__ = "Yue O Hu and Luisa W Hugerth"
__email__ = "luisa.hugerth@scilifelab.se"


def blast_parser(blastfile, maxeval, mincov, minID):
	#Save all blast hits that make the cutoff in similarity and ID
	scores = defaultdict(dict)
	with open(blastfile) as csvfile:
		reader = csv.reader(csvfile, delimiter="\t")
		for row in reader:
		##queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore
			query = row[0]
			hit = row[1]
			percid = float(row[2])
			length = int(row[3])
			hlen = float(row[9]) - float(row[8])
			evalue = float(row[10])
			cov = 100*length/hlen
			score = float(row[11])
			if (evalue <= maxeval and percid >= minID and cov >= mincov):
				scores[query][hit] = score
	return scores


def toppercent(scores, percent):
	# Keeps only the percent% high-scoring hits for each query-hit pair
	returndict = defaultdict(dict)
	#keep = scores.copy()
	numhits = len(scores.keys())
	minkeep = int(math.floor(0.01 * percent * numhits))
			#I take floor, but since lists are 0-indexed, it corresponds to ceiling
	scrs = list()
	for query, score in scores.iteritems():
		scrs.append(score)
	scrs.sort()
	if len(scrs) > minkeep:
		minscore = scrs[minkeep]
	elif len(scrs) > 0:
		minscore = scrs[-1]
	else:
		minscore = 0
	returndict = dict()
	for query, score in scores.iteritems():
		if score >= minscore:
			returndict[query] = score
	return returndict

		
def parse_taxonomy(taxfile):
	tax = defaultdict(list)
	with open(taxfile) as csvfile:
		reader = csv.reader(csvfile, delimiter="\t")
		for row in reader:
		#AB353770.1.1740_U	Eukaryota;Alveolata;Dinophyta;Dinophyceae;Dinophyceae_X;Dinophyceae_XX;Peridiniopsis;Peridiniopsis+kevei;
			tax[row[0]] = row[1]
	return tax


def lca(scores1, scores2, percent, tax):
	classdict = dict()
	for query, hit in scores1.iteritems():
		scr1 = set(hit.keys())
		scr2 = set(scores2[query].keys())
		#find the common hits of both dictionaries
		common = scr1.intersection(scr2)
		commonscores=dict()
		for goodhit in common:
			score = hit[goodhit] + scores2[query][goodhit]
			commonscores[goodhit] = score
		#get the top percent scores of this intersection
		topcommon = toppercent(commonscores, percent)
		#get the LCA for these
		classify = ''
		for query, score in topcommon.iteritems():
			if classify == '':
				classify = tax[query]
			else:
				classify = commonprefix([classify, tax[query]])
		if classify == '':
			classify = 'Unclassified;'
		#print classify
		#take longest substr ending in ;
		meaningful = re.match(".+;", classify)
		classify = meaningful.group()
		classdict[query] = classify
		#print query + "\t" + classify
	return classdict

def print_class(classify, identity):
	print "Query\tTaxonomy\tSimilarity_level"
	for query, tax in classify.iteritems():
		print query + "\t" + tax
	
def main(blast1, blast2, evalue, coverage, identity, keep, taxonomy):
##### METHOD #######
#1 filter the blast result with cutoffs 90, 97 or 99 and aligned length 250bp
#2 find the blast matches in both FWD and REV qualified items
#3 sum the score of those matches and rank them
#4 use top 5% blast annotations' Last Common Ancestor for the classification
	#parse taxonomy
	tax = parse_taxonomy(taxonomy)
	#parse blast result for each file at user cutoff
	scores1 = blast_parser(blast1, evalue, coverage, identity)
	scores2 = blast_parser(blast2, evalue, coverage, identity)
	#find the results common for forward and reverse, rank them and get the LCA of the top 5%
	classify = lca(scores1, scores2, keep, tax)
	#print out the result
	print_class(classify, identity)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Takes 2 blast output files and a taxonomy table and returns best taxonomic classification')
	parser.add_argument('-1', '--blast1', help='Blast output for forward reads')
	parser.add_argument('-2', '--blast2', help='Blast output for reverse reads')
	parser.add_argument('-e', '--evalue', nargs='?', default=1E-5, type=float, help='Maximal evalue to consider blast match. Default: %(default)f')
	parser.add_argument('-c', '--coverage', nargs='?', default=85.0, type=float, help='Minimal coverage of blast hit to consider match. Default: %(default)d per cent')
	parser.add_argument('-id', '--identity', nargs='?', default=99.0, type=float, help='Maximal residue identity of interest to consider a match. Default: %(default)d per cent')
	parser.add_argument('--keep', nargs='?', default=5.0, type=float, help='Percentage of top blast hits from each output file to check if they are paired. Default: %(default)d per cent')
	parser.add_argument('-tax', '--taxonomy', help='Annotated taxonomy for last common ancestor inference')
	args = parser.parse_args()

	main(args.blast1, args.blast2, args.evalue, args.coverage, args.identity, args.keep, args.taxonomy)

