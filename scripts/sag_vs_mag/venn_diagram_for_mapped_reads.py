# coding: utf-8

import pandas as pd
import pysam
import argparse
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context('talk')


# Group SAG reads into groups:
# 1. Reads that map to both SAG and MAG x
# 2. Reads that map to MAG but not to SAG x
# 3. Reads that map to metagenome contig longer than 1000 but not MAG and SAG x
# 4. Reads that map to metagenome contig longer than 1000 but not MAG and not to SAG x
# 5. Reads that map to metagenome contig shorter than 1000 and SAG x
# 6. Reads that map to metagenome contig shorter than 1000 but not SAG x
# 7. Reads that do not map to metagenome but to SAG x

def get_all_stats(args):
    result_l = []
    for mag_bam_file, mag, mag_contig_list, sag_bam_file, sag in zip(args.mag_bam_files, args.mag_names, args.mag_contig_lists, args.sag_bam_files, args.sag_names):
        mag_reads = set()
        long_contig_reads = set()
        metagenome_reads = set()
        date = mag.split('-')[0]

        samfile = pysam.AlignmentFile(mag_bam_file, 'rb')
        bin_contigs = pd.read_table(mag_contig_list, header=None, names=["contig_id"])
        all_mag_contigs = set(bin_contigs['contig_id'].values)
        for contig, contig_len in zip(samfile.references, samfile.lengths):
            contig_read_set = set((read.qname, read.is_read1) for read in samfile.fetch(contig))
            metagenome_reads |= contig_read_set
            if contig_len >= 1000:
                long_contig_reads |= contig_read_set
            if contig in all_mag_contigs:
                mag_reads |= contig_read_set

        sag_samfile_path = sag_bam_file
        sag_samfile = pysam.AlignmentFile(sag_samfile_path, 'rb')
        sag_reads = set((read.qname, read.is_read1) for read in sag_samfile.fetch())

        tmp_result_d = {}
        tmp_result_d['1'] = len(mag_reads & sag_reads)
        tmp_result_d['2'] = len(mag_reads - sag_reads)
        tmp_result_d['3'] = len((long_contig_reads - mag_reads) & sag_reads)
        tmp_result_d['4'] = len(long_contig_reads - mag_reads - sag_reads)
        tmp_result_d['5'] = len((metagenome_reads - long_contig_reads) & sag_reads)
        tmp_result_d['6'] = len((metagenome_reads - long_contig_reads) - sag_reads)
        tmp_result_d['7'] = len(sag_reads - metagenome_reads)
        result_l.append(tmp_result_d)

    return result_l

def plot_results(args, result_l):
    for mag, sag, sag_result_d, output_figure in zip(args.mag_names, args.sag_names, result_l, args.output_figures):

        tot_sum = 0
        for value in sag_result_d.values():
            tot_sum += value

        circle1=plt.Circle((.45,.55),.4, fill=False, color='r', lw=1.0)
        circle2=plt.Circle((0.65,0.55),.2, fill=False, color='c', lw=1.0)
        circle3=plt.Circle((0.55,0.55),.3, fill=False, color='black', lw=1.0)
        circle4=plt.Circle((.7,.3),.3, fill=False, color='b', lw=1.0)

        fig = plt.gcf()
        fig.gca().add_artist(circle1)
        fig.gca().add_artist(circle2)
        fig.gca().add_artist(circle3)
        fig.gca().add_artist(circle4)

        plt.text(.65, .45, "{:.2%}".format(sag_result_d['1'] / float(tot_sum)))
        plt.text(.65, .65, "{:.2%}".format(sag_result_d['2'] / float(tot_sum)))
        plt.text(.52, .3, "{:.2%}".format(sag_result_d['3'] / float(tot_sum)))
        plt.text(.35, .55, "{:.2%}".format(sag_result_d['4'] / float(tot_sum)))
        plt.text(.45, .2, "{:.2%}".format(sag_result_d['5'] / float(tot_sum)))
        plt.text(.25, .3, "{:.2%}".format(sag_result_d['6'] / float(tot_sum)))
        plt.text(.7, .15, "{:.2%}".format(sag_result_d['7'] / float(tot_sum)))

        plt.legend((circle1, circle2, circle3, circle4), ('Metagenome', 'MAG', 'Contigs > 1kb', 'SAG'))
        plt.title('Venn Diagram for all mapped reads. MAG: {} SAG: {}'.format(mag, sag))
        ax = plt.gca()
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid()
        ax.set_axis_bgcolor((1, 1, 1))

        fig.savefig(output_figure)

def main(args):
    if args.use_agg:
        plt.switch_backend('agg')
    result_l = get_all_stats(args)
    plot_results(args, result_l)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag_bam_files", nargs='*')
    parser.add_argument("--mag_names", nargs='*')
    parser.add_argument("--mag_contig_lists", nargs='*')
    parser.add_argument("--sag_bam_files", nargs='*')
    parser.add_argument("--sag_names", nargs='*')
    parser.add_argument("--output_figures", nargs='*')
    parser.add_argument("--use_agg", action="store_true")
    args = parser.parse_args()
    main(args)
