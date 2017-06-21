from argparse import ArgumentParser
from collections import defaultdict
import logging

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def read_bed_lines(fh):
    bed = defaultdict(list)
    for line in fh:
        line = line.rstrip()
        contig,start,end,gene = line.split("\t")
        bed[contig].append(gene)
    return bed


def read_bedfile(bedfile):
    with open(bedfile) as fh:
        bed = read_bed_lines(fh)
    return bed


def read_classifications_lines(fh,bed):
    geneassign = {}
    for line in fh:
        line = line.rstrip()
        items = line.split("\t")
        contig = items[0]
        genes = bed[contig]
        for gene in genes: geneassign[gene] = items[1:]
    return geneassign


def read_classifications(classifications,bed):
    with open(classifications) as fh:
        geneassign = read_classifications_lines(fh,bed)
    return geneassign


def write_gene_assign(geneassign):
    for gene,assign in geneassign.items():
        print(gene,'\t'.join(assign),sep="\t")


def main():
    parser = ArgumentParser()
    parser.add_argument("-b", "--bedfile", required=True,
                        help="Bed file with contig ids in first column and gene ids in 4th column")
    parser.add_argument("-c", "--classifications", required=True,
                        help="Classification output for contigs from ClassifyContigNR.py")
    args = parser.parse_args()

    logging.info("Reading bed file")
    bed = read_bedfile(args.bedfile)
    logging.info("Finished reading bed file")
    logging.info("Reading classifications file")
    geneassign = read_classifications(args.classifications,bed)
    logging.info("Finished reading classifications file")
    logging.info("Writing assignments")
    write_gene_assign(geneassign)


if __name__ == '__main__':
    main()