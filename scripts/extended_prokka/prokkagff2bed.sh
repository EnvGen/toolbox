#!/bin/bash

infile=$1

if [ "$infile" == "" ] ; then
    echo "Usage: prokkagff2bed.sh <PROKKA gff file>"
    exit 0
fi

grep -v "#" $infile | grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | awk -v OFS='\t' '{print $1, $4-1, $5,$9}'
