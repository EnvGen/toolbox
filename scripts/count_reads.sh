#!/usr/bin/env bash

#A script to count the number of reads in a gzipped fastq file and calculate average read length
# Usage ./read_count.sh <file.fastq.gz>

wc_string=`gunzip -c $1 | grep -B 1 --no-group-separator '^+$' | grep -v '^+$' | wc`
nr_chars=`echo $wc_string | cut -f 3 -d ' '`
nr_rows=`echo $wc_string | cut -f 1 -d ' '`

mean_read_length=`echo "($nr_chars - $nr_rows) / $nr_rows.0" | bc -l`

echo $mean_read_length,$nr_rows
