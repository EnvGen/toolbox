#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from collections import defaultdict
from BCBio import GFF

def read_blast_output(blastoutfile): 
    sseq_ids = []
    records = []
    with open(blastoutfile) as in_handle:
        for line in in_handle:
            line_items = line.split("\t")
            qseq = line_items[0]
            sseq = line_items[1]
            evalue = line_items[2]
            pident = line_items[3]
            send = line_items[7]
            sstart = line_items[8]
            slen = line_items[10]

            records.append({'qseqid': qseq,
                            'sseqid': sseq,
                            'evalue': float(evalue),
                            'pident': float(pident),
                            'send': float(send),
                            'sstart': float(sstart),
                            'slen': float(slen)})
            sseq_ids.append(sseq.split('|')[2])
    return records, sseq_ids

class BlastFilterer(object):
    """Class just to be able to run blast_record_ok method within a 'filter' call
    """
    def __init__(self, scovs_threshold, evalue_threshold):
        self.scovs_threshold = scovs_threshold
        self.evalue_threshold = evalue_threshold

    def blast_record_ok(self, blast_record):
        """
        Checks for one blast record if it passes threshold filters
        """
        evalue_above_threshold = blast_record['evalue'] >= self.evalue_threshold

        alignment_length_in_subject = abs(blast_record['send'] - blast_record['sstart']) + 1
        percent_seq_covered = (alignment_length_in_subject / blast_record['slen']) * 100.0
        seq_cov_above_threshold =  percent_seq_covered >= self.scovs_threshold
    
        return evalue_above_threshold and seq_cov_above_threshold

class GFF_feature_rec(object):
    """
    Class to represent a feature record from a row in a gff file
    """
    def __init__(self, row):
        entries = row.split('\t')
        self.seqid = entries[0]
        self.source = entries[1]
        self.type = entries[2]
        self.start = entries[3]
        self.end = entries[4]
        self.score = entries[5]
        self.strand = entries[6]
        self.phase = entries[7]
        self.attributes = entries[8]
        
        attributes = self.attributes.split(';')
        ID = attributes[0]
        if ID.startswith('ID='):
            self.ID = ID.split('=')[1]
        

def extend_record(row, feature_blast_hit, prokka=True):
    """
    extends record rec with the feature feature_blast_hit if 
    not already present in record.
    """
    rec = GFF_feature_rec(row)
    if rec.ID:
        new_xref = None
        if prokka and rec.ID in feature_blast_hit:
            new_xref = feature_blast_hit[rec.ID]
        elif not prokka:
            prodigal_rec_id = rec.seqid + '_' + rec.ID.split('_')[-1] 
            if prodigal_rec_id in feature_blast_hit:
                new_xref = feature_blast_hit[prodigal_rec_id]
        if new_xref is not None:
            row += ";Dbxref=" + new_xref
    return row

def translate_from_cdd(feature_blast_hits, cdd_all_file, include_evalue, include_pident):
    with open(cdd_all_file, 'r') as cf:
         cddid_d = dict([(row.split('\t')[0], row.split('\t')[1].strip()) for row in cf.readlines()])
         preformatted_dict = {}
         for k, v in feature_blast_hits.iteritems():
             new_val = cddid_d[v['sseqid'].split('|')[-1]]
             if include_evalue:
                 new_val += ",evalue:{0}".format(v['evalue'])
             if include_pident:
                 new_val += ",pident:{0}".format(v['pident'])
             preformatted_dict[k] = new_val
    return preformatted_dict 

def main(blastoutfile, gff_file, scovs_threshold, evalue_threshold, cddid_all_file, include_evalue, include_pident, not_prokka):
    # Read the blast output
    blast_records, sseq_ids = read_blast_output(blastoutfile)

    # Filter on evalue and scovs-threshold
    blast_filter = BlastFilterer(scovs_threshold, evalue_threshold)
    blast_records = filter(blast_filter.blast_record_ok, blast_records)
    
    # Convert blast records to dict indexed on qseqid
    feature_blast_hits =  dict([(blast_r['qseqid'], blast_r) for blast_r in blast_records])    

    # We only want the translated cdd in the resulting genbank file
    translated_blast_hits = translate_from_cdd(feature_blast_hits, cddid_all_file, include_evalue, include_pident)

    new_gff = []
    include_rest = False
    with open(gff_file, 'r') as gff_fh:
        for i, row in enumerate(gff_fh):
            row = row.strip()
            if include_rest:
                new_gff.append(row)
            elif row.startswith('#'): #The initial header lines
                new_gff.append(row)
                
            elif row.startswith('>'): #The ending fasta seqs
                new_gff.append(row)
                include_rest = True
            else:
                new_gff.append(extend_record(row, translated_blast_hits, prokka=not(not_prokka)))

    sys.stdout.write('\n'.join(new_gff))
     


if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument('-b', '--blastoutfile', required=True,
           help=('Output of rpsblast run, assumed to be in tabular format whith '
               'columns: qseqid sseqid evalue pident score qstart qend sstart send length slen. '
               'The contigs ids are assumed to be recoverable by removing the last underscore '
               'and the characters following it from the qseqid column.' ))
   parser.add_argument('-g', '--gfffile', required=True,
           help=('General Feature Format file to be extended on stdout.'))
   parser.add_argument('-s', '--scovs-threshold', type=float, default=50.0,
           help='Threshold covered in percent, default=50.0')
   parser.add_argument('-e', '--evalue-threshold', type=float, default=0.0,
           help='Threshold evalue, default=0.0')
   parser.add_argument('--cddid_all_file', required=True,
           help = ('A table listing all entries in CDD together with its '
           'accessions. Used to find e.g. Pfam and TIGRFAM ids given the '
           'the CDD ids. '
           'Downloaded from: ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz'))
   parser.add_argument('--include-evalue', action='store_true',
           help = ('Use this tag if evalue is supposed to be in the db_xref '
               'tag along side the translated CDD id.'))
   parser.add_argument('--include-pident', action='store_true',
           help = ('Use this tag if pident is supposed to be in the db_xref ' 
               'tag along side the translated CDD id.'))
   parser.add_argument('--not_prokka', action='store_true',
           help = ('Use this tag if the gff file and the proteins are not '
               'generated by prokka.'))
   args = parser.parse_args()

   main(args.blastoutfile, 
        args.gfffile, 
        args.scovs_threshold, 
        args.evalue_threshold, 
        args.cddid_all_file,
        args.include_evalue,
        args.include_pident,
        args.not_prokka)

