# /usr/bin/env python3
import os
import argparse
import sys
import io
import gzip

# Author: Luis Fernando Delgado

usage = 'usage: create_gff.py [-i] [-o]'
description = 'This program creates a "fake" gff file from representative genes'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', metavar="input fasta file", dest='i', help='nucleotide fasta. e.g, .fna or .fna.gz (interleaved or not)', required=True)
parser.add_argument('-o', metavar="output file name", dest='o', help='file name output, .gff', required=True)
args = parser.parse_args()

'''
opening sequences:
>co_assembly_k127_1_1 # 2 # 208 # 1 # ID=1_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.609
AGCCATTGCTGCTGGCAGATGCGACATTGCGCTGGCAGACGAGCCATCTCGCTGACAGACGAGCCATTGC
TGCTGGCAGACGAGCCATTGCTGCTGGCAGATGCGACATTGCGCTGGCAGACGAGCCATCTCGCTGACAG
ACGAGCCATTGCTGCTGGCAGACGAGCCATTGCTGCTGGCAGATGCGACATTGCGCTGGCAGACGAG

to be converted to:
co_assembly_k127_1_1     Prodigal_v2.6.3 CDS     1   207 . + . ID=1_1@co_assembly_k127_1_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.609
'''
# Open input fasta file
if args.i.endswith(".gz"):
    input_file = gzip.open(args.i, 'rb')
    fin = io.TextIOWrapper(input_file, encoding='utf-8')
else:
    fin = open(args.i, "r")

length = 0

with open(args.o, "w") as fout:
    for line in fin:
        line = line.rstrip()
        if line.startswith(">"):   # fasta header
            if length != 0:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(name, "Prodigal_v2.6.3", "CDS", "1", length, ".", direct, ".", rest), file=fout)
                length = 0

            head = line.split(" ")
            name = head[0][1:]    # co_assembly_k127_1_1
            direction = head[6]   # 1

            if direction.startswith("-"):
                direct = "-"
            else:
                direct = "+"   # +

            fp = head[8]     # ID=1_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.609
            part1 = fp.split(";")[0]   # ID=1_1
            part2 = fp.split(";")[1:]   # partial=11 start_type=Edge rbs_motif=None rbs_spacer=None gc_cont=0.609
            ID = part1+"@"+name   # ID=1_1@co_assembly_k127_1_1
            rest = str(ID)+";"+str(";".join(part2))   # ID=1_1@co_assembly_k127_1_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.609

        else:   # Sequence
            length += len(line)
# printing out the last sequence
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(name, "Prodigal_v2.6.3", "CDS", "1", length, ".", direct, ".", rest), file=fout)

# closing input file
if args.i.endswith(".gz"):
    input_file.close()
fin.close()
