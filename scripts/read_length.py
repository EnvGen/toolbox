#!/usr/bin/env python
import argparse

def main(args):
    total_length = 0
    with open(args.input_fastq) as ifh:
        for i, line in enumerate(ifh):
            if i % 4 == 1:
                line.strip()
                total_length+=len(line)

    if args.sample_name:
        sample_name = args.sample_name
    else:
        sample_name = args.input_fastq
    print("{}\t{}".format(sample_name, total_length / float((i+1)/4)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fastq")
    parser.add_argument("--sample_name")

    args = parser.parse_args()
    main(args)
