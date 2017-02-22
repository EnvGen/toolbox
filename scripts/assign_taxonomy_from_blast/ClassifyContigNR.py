import argparse
import logging
import operator
import os
from collections import Counter
from collections import defaultdict

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

# These are identities normalized with query coverage:
MIN_IDENTITY_TAXA = (0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

def get_id(subjectId, accession_mode):
    '''Returns the identifier for the query'''
    if accession_mode:
        return subjectId
    else:
        return subjectId.split("|")[1]


def calculate_fHit(percIdentity, queryEnd, queryStart, qLength):
    alnLength_in_query = abs(int(queryEnd) - int(queryStart)) + 1
    fHit = float(alnLength_in_query) / qLength
    fHit *= float(percIdentity) / 100.0
    fHit = min(1.0, fHit)
    return fHit


def test_fHit_calculation():
    percIdentity = 90.0
    queryEnd = 309
    queryStart = 10
    qLength = 600
    ## 309-10+1 = 300, 300/600 = 0.5, 0.5*(90/100) = 0.9*0.5 = 0.45
    assert calculate_fHit(percIdentity, queryEnd, queryStart, qLength) == 0.45


def read_blast_lines(fh, lengths, accession_mode, min_id):
    # k191_83_2       gi|973180054|gb|KUL19018.1|     71.2    73      21      0       9       81      337     409     6.6e-24 118.2
    # queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore

    matches = defaultdict(list)
    gids = Counter()
    for line in fh:
        line = line.rstrip()
        (queryId, subjectId, percIdentity, alnLength, mismatchCount,
         gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore) = line.split("\t")

        ## Get the accession or GenInfo Identifier (gi)
        gid = get_id(subjectId, accession_mode)
        ## Get length of query
        qLength = lengths[queryId]
        ## Calculate percent identity normalized to alignment length
        fHit = calculate_fHit(percIdentity, queryEnd, queryStart, qLength)

        if float(percIdentity) >= min_id:
            matches[queryId].append((gid, fHit))
            gids[gid] += 1
    return (matches, gids)


def read_blast_input(blastinputfile, lengths, min_id, accession_mode=False):
    with open(blastinputfile) as fh: (matches, gids) = read_blast_lines(fh, lengths, accession_mode, min_id)
    return (matches, gids.keys())


def run_map_back(mapBack, tokens):
    for depth in range(6, 0, -1):
        if tokens[depth] not in mapBack[depth] and tokens[depth] != 'None':
            for depth2 in range(depth - 1, -1, -1):
                mapBack[depth][tokens[depth]].append(tokens[depth2])
    return mapBack


def read_lineage_lines(fh):
    mapping = {}
    mapBack = defaultdict(lambda: defaultdict(list))
    for line in fh:
        line = line.rstrip()
        tokens = line.split("\t")
        (taxaid, domainid, phylumid, classid, orderid, familyid, genusid, speciesid) = tokens
        ## Map taxaid to the lineage
        mapping[int(taxaid)] = [domainid, phylumid, classid, orderid, familyid, genusid, speciesid];
        tokens.pop(0)  ## Remove taxaid from tokens
        ## Iterate over the ranks and map each rank to the full previous lineage
        mapBack = run_map_back(mapBack, tokens)
    return (mapping, mapBack)


def read_lineage_file(lineage_file):
    with open(lineage_file) as fh: (mapping, mapBack) = read_lineage_lines(fh)
    return (mapping, mapBack)


def read_query_length_lines(fh):
    lengths = {}
    for line in fh:
        line = line.rstrip()
        (queryid, length) = line.split("\t")
        lengths[queryid] = float(length)
    return lengths


def read_query_length_file(query_length_file):
    with open(query_length_file) as fh: lengths = read_query_length_lines(fh)
    return lengths


def map_gids_binary(gids, mapping_file):
    t = open(mapping_file, 'r')
    c = 0
    size = long(os.path.getsize(mapping_file))
    mapping = {}

    for gis in gids:  # binary search #########################################
        c = c + 1
        # print gis
        found = False
        offset = 0
        chunk = size
        pos = chunk / 2
        while found == False and chunk > 0:
            chunk = chunk / 2
            t.seek(pos)
            t.readline()
            entry = t.readline().split("\t")
            if entry[0]:
                filegi = int(entry[0])
                filetax = int(entry[1].rstrip("\n"))

                if filegi == int(gis):
                    answer = filetax
                    found = True
                elif filegi > int(gis):
                    pos = offset + (chunk / 2)
                elif filegi < int(gis):
                    offset = offset + chunk
                    pos = pos + (chunk / 2)

        if found == False:
            answer = -1

        mapping[gis] = answer

    return mapping


def map_accessions(accs, fh):
    mappings = dict([(acc, -1) for acc in accs])
    for i, line in enumerate(fh):
        if i == 0: continue
        _, acc_ver, taxid, _ = line.split("\t")
        # Only add taxids for the given acc
        if acc_ver in mappings:
            mappings[acc_ver] = int(taxid)
    return mappings


def read_accessions_file(accs, mapping_file):
    with open(mapping_file) as fh: mappings = map_accessions(accs, fh)
    return mappings


def calculate_gene_length(start, end):
    return abs(int(end) - int(start)) + 1


def test_calculate_gene_length():
    assert calculate_gene_length(10, 1) == 10
    assert calculate_gene_length(1, 10) == 10


def read_bed_lines(fh):
    contigGenes = defaultdict(list)
    contigLengths = Counter()
    lengths = {}
    for line in fh:
        line = line.rstrip()
        contig, start, end, gene = line.split("\t")
        contigGenes[contig].append(gene)
        length = calculate_gene_length(start, end)
        lengths[gene] = length
        contigLengths[contig] += length
    return (lengths, contigGenes, contigLengths)


def read_bed_file(bedfile):
    with open(bedfile) as fh: (lengths, contigGenes, contigLengths) = read_bed_lines(fh)
    return lengths, contigGenes, contigLengths


def calculate_taxa_weight(fHit, min_id_taxa):
    ## If min_id_taxa = 0.95, fHit needs to be > 0.95 for this depth
    weight = (fHit - min_id_taxa) / (1.0 - min_id_taxa)
    weight = max(weight, 0.0)
    return weight


def test_calculate_taxa_weight():
    fHit = 0.75
    min_id_taxa = 0.5
    assert calculate_taxa_weight(fHit, min_id_taxa) == 0.5


def collate_gene_hits(matchs, mapping, lineages):
    collate_hits = list()
    for depth in range(7): collate_hits.append(Counter())

    added_matches = set()
    ## Iterate the hits (protein accessions) and normalized identity, sorted by decreasing normalized percent identity
    for (match, fHit) in sorted(matchs, key=lambda x: x[1], reverse=True):
        if mapping[match] <= -1: continue
        ## Get the taxid and lineage for the protein accession
        tax_id = mapping[match]

        ## Only add best hit per species
        if tax_id in added_matches: continue
        added_matches.add(tax_id)

        if tax_id not in lineages:
            logging.warning("Taxa id {} is missing from lineage file".format(tax_id))
            continue

        hits = lineages[tax_id]
        for depth in range(7):
            if hits[depth] != "None":
                ## Calculate the normalized weight for each depth
                weight = calculate_taxa_weight(fHit, MIN_IDENTITY_TAXA[depth])
                if weight > 0.0: collate_hits[depth][hits[depth]] += weight  # could put a transform in here
    return collate_hits


def make_assignment(geneAssign, gene, collate_hits, mapBack, min_fraction):
    for depth in range(6, -1, -1):
        collate = collate_hits[depth]
        dWeight = sum(collate.values())
        sortCollate = sorted(collate.items(), key=operator.itemgetter(1), reverse=True)
        nL = len(collate)  ## Number of hits collated at depth
        if nL > 0:
            dP = 0.0
            ## dWeight is the sum of weights at this depth
            if dWeight > 0.0:
                ## dP is the weight of the best hit at depth, divided by the total weight at depth
                dP = float(sortCollate[0][1]) / dWeight
                ## If the dP for the best hit is larger than min_fraction then assign at this depth and map back
                if dP > min_fraction:
                    geneAssign[gene][depth] = (sortCollate[0][0], dP)
                    ## Get back lineage for the best assignment
                    assignBack = mapBack[depth][sortCollate[0][0]]
                    depth2 = depth - 1
                    for assignB in assignBack:
                        geneAssign[gene][depth2] = (assignB, 1.0)
                        depth2 -= 1
                    break

                else:
                    geneAssign[gene][depth] = ('Unclassified', -1.0)
            ## If dWeight is 0.0, set 'Unclassified'
            else:
                geneAssign[gene][depth] = ('Unclassified', -1.0)
        ## If no hits collated at depth, set 'Unclassified'
        else:
            geneAssign[gene][depth] = ('Unclassified', -1.0)
    return geneAssign


def assign_unclassified():
    d = {}
    for depth in range(7): d[depth] = ('Unclassified',-1.0)
    return d

def assign_taxonomy(matches, mapping, mapBack, lineages, lengths, contigGenes, contigLengths, min_fraction):
    geneAssign = defaultdict(dict)
    contigAssign = defaultdict(dict)
    for contig, genes in contigGenes.items():
        collate_hits = list()
        for depth in range(7):
            collate_hits.append(Counter())
        for gene in genes:
            ## Perform taxonomic assignments for all genes on contig
            try:
                matchs = matches[gene]
                collated_gene_hits = collate_gene_hits(matchs, mapping, lineages)
                geneAssign = make_assignment(geneAssign, gene, collated_gene_hits, mapBack, min_fraction)
            except KeyError:
                geneAssign[gene] = assign_unclassified()
                continue
            ## Then add assignments to the contig for collating
            for depth in range(7):
                try:
                    (assignhit, genef) = geneAssign[gene][depth]
                except KeyError:
                    continue
                if assignhit != 'Unclassified':
                    collate_hits[depth][assignhit] += lengths[gene]

        # Contigs are assigned a taxonomy based on LCA for
        # all genes assigned on at least kingdom level.
        dWeight = sum(collate_hits[0].values())
        for depth in range(7):
            collate = collate_hits[depth]
            sortCollate = sorted(collate.items(), key=operator.itemgetter(1), reverse=True)
            nL = len(collate)
            if nL > 0:
                dP = 0.0
                if dWeight > 0.0:
                    dP = float(sortCollate[0][1]) / dWeight
                    if dP > min_fraction:
                        contigAssign[contig][depth] = (sortCollate[0][0], dP, sortCollate[0][1])
                    else:
                        contigAssign[contig][depth] = ('Unclassified', 0., 0.)
                else:
                    contigAssign[contig][depth] = ('Unclassified', 0., 0.)
            else:
                contigAssign[contig][depth] = ('Unclassified', 0., 0.)

    return contigAssign,geneAssign


def write_gene_assigns(output_dir, geneAssign):
    with open(output_dir + "_genes.csv", "w") as assign_file, open(output_dir + "_genes.supports.csv",
                                                                   "w") as support_file:
        for gene in geneAssign.keys():
            assign_file.write('%s' % gene)
            support_file.write('%s' % gene)
            for depth in range(7):
                (assign, p) = geneAssign[gene][depth]
                assign_file.write(',%s' % assign)
                support_file.write(',%.3f' % p)
            assign_file.write('\n')
            support_file.write('\n')
            assign_file.flush()
            support_file.flush()


def write_contig_assigns(output_dir, contigAssign, contigLengths):
    with open(output_dir + "_contigs.csv", "w") as assign_file, open(output_dir+"_contigs.supports.csv", "w") as support_file:
        for contig in contigAssign.keys():
            assign_file.write('%s,%f' % (contig, contigLengths[contig]))
            for depth in range(7):
                (assign, p, dF) = contigAssign[contig][depth]
                dFN = dF / contigLengths[contig]
                assign_file.write(',%s'%assign)
                support_file.write(',%.3f/%.3f' % (p, dFN))
            assign_file.write('\n')
            support_file.write('\n')
            assign_file.flush()
            support_file.flush()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_input_file",
                        help="directory with blast 6 matches to taxaid database *.b6")
    # parser.add_argument("query_length_file", help="tab delimited file of query lengths")
    parser.add_argument('-g', '--gid_taxaid_mapping_file',
                        help="mapping from gid to taxaid gzipped")
    parser.add_argument('-a', '--acc_taxaid_mapping_file',
                        help="mapping from accession to taxaid gzipped")
    parser.add_argument('-b', '--bed_file', required=True,
                        help="bed file defining genes on contigs")
    parser.add_argument('-l', '--lineage_file', required=True,
                        help="text taxaid to lineage mapping")
    parser.add_argument('-f', '--min_fraction', default=0.5, type=float,
                        help="Minimum fraction of weights needed to assign to a particular level, default 0.5. 0.9 would be more strict.")
    parser.add_argument('-i', '--min_id', default=40.0, type=float,
                        help="Minimum allowed percent identity to parse a hit")
    parser.add_argument('-o', '--output_dir', type=str, default="output",
                        help=("string specifying output directory and file stubs"))

    args = parser.parse_args()

    if args.gid_taxaid_mapping_file and args.acc_taxaid_mapping_file:
        raise Exception(
            "Both gid_taxaid_mapping_file and acc_taxaid_mapping_file are given, but only one at a time is allowed")
    elif args.gid_taxaid_mapping_file:
        accession_mode = False
    else:
        accession_mode = True

    (lengths, contigGenes, contigLengths) = read_bed_file(args.bed_file)
    logging.info("Finished reading bed file")

    (matches, gids) = read_blast_input(args.blast_input_file, lengths, args.min_id, accession_mode)
    logging.info("Finished reading in blast results file")

    (lineages, mapBack) = read_lineage_file(args.lineage_file)
    logging.info("Finished reading in lineage file")

    if accession_mode:
        mapping = read_accessions_file(gids, args.acc_taxaid_mapping_file)
    else:
        mapping = map_gids_binary(gids, args.gid_taxaid_mapping_file)
    logging.info("Finished loading taxaid map file")

    logging.info("Assigning taxonomy")
    contigAssign, geneAssign = assign_taxonomy(matches, mapping, mapBack, lineages, lengths, contigGenes, contigLengths, args.min_fraction)
    write_gene_assigns(args.output_dir, geneAssign)
    write_contig_assigns(args.output_dir, contigAssign, contigLengths)
    #geneAssign = assign_taxonomy_to_genes(matches, mapping, mapBack, lineages, args.min_fraction)
    #logging.info("Finished assigning taxonomy to genes")
    #write_gene_assigns(args.output_dir, geneAssign)

    #contigAssign = assign_taxonomy_to_contigs(geneAssign, lengths, contigGenes, contigLengths, args.min_fraction)
    #logging.info("Finished assigning taxonomy to contigs")
    #write_contig_assigns(args.output_dir, contigAssign, contigLengths)


if __name__ == "__main__":
    main()
