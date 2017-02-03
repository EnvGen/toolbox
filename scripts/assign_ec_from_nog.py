#!/usr/bin/env python
import argparse
from collections import defaultdict

def main(args):
    nog_to_ec = defaultdict(set)
    with open(args.nog_to_ec, 'r') as ifh:
        for line in ifh:
            line = line.strip()
            fields = line.split('\t')
            try:
                nog_to_ec[fields[0]].add(fields[1])
            except:
                print(fields)
                raise

    with open(args.NOG_annotation, 'r') as ifh:
        for line in ifh:
            line = line.strip()
            fields = line.split('\t')
            if fields[1] in nog_to_ec:
                for ec in nog_to_ec[fields[1]]:
                    print("\t".join([fields[0], ec] + fields[2:]))
	
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("NOG_annotation")
    parser.add_argument("nog_to_ec")
    args = parser.parse_args()
    main(args)
