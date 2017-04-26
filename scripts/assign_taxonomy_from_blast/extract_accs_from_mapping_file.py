#!/usr/bin/env python

from argparse import ArgumentParser
import sys

def read_accession_lines(fh):
    accs = {}
    for line in fh:
        line = line.rstrip()
        accs[line] = ""
    return accs


def read_accessions_file(f):
    with open(f, 'r') as fh: accs = read_accession_lines(fh)
    return accs


def read_map_lines(fh,accs):
    total = len(accs)
    found = 0
    mapping = {}
    not_found = set()
    for i,line in enumerate(fh):
        if i==0: continue
        line = line.rstrip()
        acc,acc_ver,taxid,gi = line.split("\t")
        if acc_ver in accs:
            mapping[acc_ver] = [acc,acc_ver,taxid,gi]
            found+=1
            if found == total: break
    if found < total:
        not_found = set(accs.keys()).difference(set(mapping.keys()))
    return (mapping,not_found)


def extract_mappings(accs,mapfile):
    with open(mapfile, 'r') as fh: (mapping, not_found) = read_map_lines(fh,accs)
    return (mapping, not_found)


def write_mapping(mapping):
    print("accession","accession.version","taxid","gi",sep="\t")
    for key,values in mapping.items():
        print(values[0],values[1],values[2],values[3],sep="\t")


def main():
    parser = ArgumentParser()
    parser.add_argument("accession_file", type=str,
                        help="File with accession numbers to limit mapping file to")
    parser.add_argument("mapping_file", type=str,
                        help="Protein accession to taxid mapping file")

    args = parser.parse_args()

    accs = read_accessions_file(args.accession_file)

    (mapping, not_found) = extract_mappings(accs,args.mapping_file)

    write_mapping(mapping)

    if len(not_found) > 0:
        sys.stderr.write("WARNING: Taxid not found for "+str(len(not_found))+" accessions\n")
        for acc in not_found: sys.stderr.write(acc+"\n")


if __name__ == '__main__':
    main()
