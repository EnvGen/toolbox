Extended Prokka
===============
Extended Prokka is a workflow where a modified version of the [Prokka pipeline](https://github.com/tseemann/prokka) is used first and then extended with additional annotations. The [modified version](https://github.com/EnvGen/prokka) allows genes to run off edges.

Scripts
-------

###fake_tbl2asn.py
A dummy script, only producing empty files. This script might be useful if you don't need the genbank file ('.gbk') and you have a very large input file. Using a large input file we experienced that a majority of the execution time went into the execution of tbl2asn and since we were only interested in the gff file, we created this script. To use it, rename it to tbl2asn and make sure its found before the real tbl2asn in your $PATH.

###extend_gff.py
A script to extend a gff file with all annotations found in a blast output file. We use this to for example always get all COG annotations for an input file instead of the default behaviour of prokka where it for example only reports the UniProt hit if a contig has one of those.

###get.prokka.annotation.sh
A script to grep after any special annotation in the output files of extended prokka. The output is a tab separated list of which genes are annotated as which group.

###collate.annotations.py
A script to summarize all output files from get.prokka.annotation.sh into one table per annotation type.
