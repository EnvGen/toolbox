#!/usr/bin/env python
"""
Parses a large number of dnadiff reports to only report the most significant hits

Only hits with more than 10% of the query bases aligned is reported.

output_folder/fastaname1_vs_fastaname2/
output_folder/fastaname1_vs_fastaname3/

etc

"""
import argparse
import re
import os
from os.path import join as ospj

class MUMmerReport(object):
    """Represents .report file from MUMmer's dnadiff. Stores TotalBases and
    AlignedBases (%) stats."""
    def __init__(self, report_file):
        self.report_file = report_file

        with open(report_file) as f:
            one_to_one_parsed = False

            for line in f:
                if line.startswith("TotalBases"):
                    self.tot_bases = [int(b) for b in line.split()[1:]]
                if line.startswith("AlignedBases"):
                    # store percentage
                    self.aligned_perc = [float(p) for p in
                            re.findall(r'\((.*?)\%\)', line)]
                    self.aligned_bases = [int(n) for n in 
                            re.findall(r'\W([0-9]*)\(', line)]
                if not one_to_one_parsed and line.startswith("AvgIdentity"):
                    self.avg_identity = [float(p) for p in line.split()[1:]][0]
                    one_to_one_parsed = True
        self.aligned_perc_ref = self.aligned_perc[0]
        self.aligned_perc_mag = self.aligned_perc[1]
        self.aligned_bases_ref = self.aligned_bases[0]
        self.aligned_bases_mag = self.aligned_bases[1]

def main(args):
    # Get fasta names

    # Get basename from fasta files and see if those are unique
    fasta_names_ref = [".".join(os.path.basename(f).split(".")[0:-1]) for f in
            args.fasta_files_ref]
    fasta_names_mag = [os.path.basename(f).split(".")[0] for f in
            args.fasta_files_mag]

    print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format("ref_fasta_name", "mag_fasta_name", "aligned_bases_ref", "aligned_perc_ref", "aligned_bases_mag", "aligned_perc_mag", "avg_identity"))
    for ref_fasta_name in fasta_names_ref:
        for mag_fasta_name in fasta_names_mag:
            repfile = ospj(args.input_dir, "{fn1}_vs_{fn2}.report".format(
                fn1=ref_fasta_name, fn2=mag_fasta_name))
            mumr = MUMmerReport(repfile)
            if mumr.aligned_perc_mag >= args.min_coverage:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(ref_fasta_name, mag_fasta_name, mumr.aligned_bases_ref, mumr.aligned_perc_ref, mumr.aligned_bases_mag, mumr.aligned_perc_mag, mumr.avg_identity))

if __name__ == "__main__":
    """Return input arguments using argparse"""
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_dir")
    parser.add_argument("--fasta_files_ref", nargs='+')
    parser.add_argument("--fasta_files_mag", nargs='+')
    parser.add_argument("--min_coverage", type=float, default=10)
    args = parser.parse_args()
    main(args)
