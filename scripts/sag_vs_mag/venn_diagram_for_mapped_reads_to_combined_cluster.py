# coding: utf-8

import numpy as np
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
# 8. Reads that do not map to metagenome neither to SAG

def get_all_stats(args):
    result_l = []

    max_len = 0
    if args.sag_new_bam:
        arg_lists = [args.mag_bam_files, args.mag_names, args.mag_contig_lists, args.sag_bam_files, args.sag_names, args.sag_unmapped_bam_files, args.sag_new_bam]
    else:
        arg_lists = [args.mag_bam_files, args.mag_names, args.mag_contig_lists, args.sag_bam_files, args.sag_names, args.sag_unmapped_bam_files]

    for arg_list in arg_lists:
        if len(arg_list) > max_len:
            max_len = len(arg_list)
    new_arg_lists = []
    for arg_list in arg_lists:
        if len(arg_list) == max_len:
            new_arg_lists.append(arg_list)
        else:
            new_arg_lists.append(max_len*arg_list)

    for arg_list in zip(*new_arg_lists):
        if args.sag_new_bam:
            mag_bam_file, mag, mag_contig_list, sag_bam_file, sag, sag_unmapped_bam_file, sag_new_bam = arg_list
        else:
            mag_bam_file, mag, mag_contig_list, sag_bam_file, sag, sag_unmapped_bam_file = arg_list
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
        sag_mapped_reads = set((read.qname, read.is_read1) for read in sag_samfile.fetch())

        if args.sag_new_bam:
            with pysam.AlignmentFile(sag_new_bam, "wb", header=header) as outf:
                for read in sag_samfile.fetch():
                    if (read.qname, read.is_read1) in metagenome_reads:
                        if (read.qname, read.is_read1) in mag_reads:
                            read.tags = (("RG", "MAG_mapped"))
                        elif (read.qname, read.is_read1) in long_contig_reads:
                            read.tags = (("RG", "Long_non_MAG_contig_mapped"))
                        else:
                            read.tags = (("RG", "Short_metagenome_contig_mapped"))
                    else:
                        read.tags = (("RG", "Metagenome_unmapped"))
                    outf.write(read)

        sag_unmapped_samfile = pysam.AlignmentFile(sag_unmapped_bam_file, 'rb')
        sag_unmapped_reads = set((read.qname, read.is_read1) for read in sag_unmapped_samfile)

        sag_all_reads = sag_mapped_reads | sag_unmapped_reads
        # Use intersection with the sag_mapped_reads for all sets
        # this will remove duplicates as have been done 
        # against the sag sequences.
        tmp_result_d = {}
        tmp_result_d['1'] = len(mag_reads & sag_mapped_reads)
        tmp_result_d['2'] = len((mag_reads - sag_mapped_reads) & sag_all_reads)
        tmp_result_d['3'] = len((long_contig_reads - mag_reads) & sag_mapped_reads) # sag_mapped_reads is a subset of sag_all_reads
        tmp_result_d['4'] = len((long_contig_reads - mag_reads - sag_mapped_reads) & sag_all_reads)
        tmp_result_d['5'] = len((metagenome_reads - long_contig_reads) & sag_mapped_reads) # sag_mapped_reads is a subset of sag_all_reads
        tmp_result_d['6'] = len(((metagenome_reads - long_contig_reads) - sag_mapped_reads) & sag_all_reads)
        tmp_result_d['7'] = len(sag_mapped_reads - metagenome_reads) # sag_mapped_reads is a subset of sag_all_reads
        tmp_result_d['8'] = len(sag_unmapped_reads - metagenome_reads) # sag_unmapped_reads is a subset of sag_all_reads

        result_l.append(tmp_result_d)


    if args.table_output:
        t_o_df = pd.DataFrame(result_l, index=zip(new_arg_lists[1], new_arg_lists[4]))
        t_o_df.to_csv(args.table_output, sep='\t')

    if args.summarize_in_one_plot:
        summary_result = {}
        summary_std_dev = {}
        for stat_nr in range(1,9):
            stat_s = str(stat_nr)
            data_a = np.array([tmp_result_d[stat_s] for tmp_result_d in result_l])
            summary_result[stat_s] = data_a.mean()
            summary_std_dev[stat_s] = data_a.std()
        result_l = [summary_result]
        return result_l, None

    return result_l, None

def plot_results(title, sag_result_d, output_figure, std_dev):
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

    if std_dev is None:
        plt.text(.65, .45, "{:.2%}".format(sag_result_d['1'] / float(tot_sum)))
        plt.text(.65, .65, "{:.2%}".format(sag_result_d['2'] / float(tot_sum)))
        plt.text(.52, .3, "{:.2%}".format(sag_result_d['3'] / float(tot_sum)))
        plt.text(.35, .55, "{:.2%}".format(sag_result_d['4'] / float(tot_sum)))
        plt.text(.45, .2, "{:.2%}".format(sag_result_d['5'] / float(tot_sum)))
        plt.text(.25, .3, "{:.2%}".format(sag_result_d['6'] / float(tot_sum)))
        plt.text(.7, .15, "{:.2%}".format(sag_result_d['7'] / float(tot_sum)))
        plt.text(.1, .1, "{:.2%}".format(sag_result_d['8'] / float(tot_sum)))
    else:
        plt.text(.65, .45, "{:.2%} ({:.2})".format(sag_result_d['1'] / float(tot_sum), float(std_dev['1'])))
        plt.text(.65, .65, "{:.2%} ({:.2})".format(sag_result_d['2'] / float(tot_sum), std_dev['2']))
        plt.text(.52, .3, "{:.2%} ({:.2})".format(sag_result_d['3'] / float(tot_sum), std_dev['3']))
        plt.text(.35, .55, "{:.2%} ({:.2})".format(sag_result_d['4'] / float(tot_sum), std_dev['4']))
        plt.text(.45, .2, "{:.2%} ({:.2})".format(sag_result_d['5'] / float(tot_sum), std_dev['5']))
        plt.text(.25, .3, "{:.2%} ({:.2})".format(sag_result_d['6'] / float(tot_sum), std_dev['6']))
        plt.text(.7, .15, "{:.2%} ({:.2})".format(sag_result_d['7'] / float(tot_sum), std_dev['7']))
        plt.text(.1, .1, "{:.2%} ({:.2})".format(sag_result_d['8'] / float(tot_sum), std_dev['8']))
    plt.legend((circle1, circle2, circle3, circle4), ('Metagenome', 'MAG', 'Contigs > 1kb', 'SAG'))
    plt.title(title)
    ax = plt.gca()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid()
    ax.set_axis_bgcolor((1, 1, 1))

    fig.savefig(output_figure)

def plot_all_results(args, result_l, std_dev):
    if args.summarize_in_one_plot:
        title = 'Venn Diagram for all mapped reads. {}'.format(args.name)
        plot_results(title, result_l[0], args.output_figures[0], std_dev)
    else:
        for mag, sag, sag_result_d, output_figure in zip(args.mag_names, args.sag_names, result_l, args.output_figures):
            title = 'Venn Diagram for all mapped reads. MAG: {} SAG: {}'.format(mag, sag)
            plot_results(title, sag_result_d, output_figure, None)

def main(args):
    if args.use_agg:
        plt.switch_backend('agg')
    result_l, std_dev = get_all_stats(args)
    plot_all_results(args, result_l, std_dev)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag_bam_files", nargs='*')
    parser.add_argument("--mag_names", nargs='*')
    parser.add_argument("--mag_contig_lists", nargs='*')
    parser.add_argument("--sag_bam_files", nargs='*')
    parser.add_argument("--sag_names", nargs='*')
    parser.add_argument("--sag_unmapped_bam_files", nargs='*')
    parser.add_argument("--sag_new_bam", nargs='*')
    parser.add_argument("--output_figures", nargs='*')
    parser.add_argument("--use_agg", action="store_true")
    parser.add_argument("--summarize_in_one_plot", action="store_true", help="Use this tag if all input files will generate one plot")
    parser.add_argument("--name", help="A string to be used as part of the title for summarized plots")
    parser.add_argument("--table_output", help="File to which a tsv representation of the data will be written") 
    args = parser.parse_args()
    main(args)
