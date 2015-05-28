#!/bin/bash

infile="$1"
type="$2"

if [ "$infile" == "" ] ; then
    echo -e "Usage: get.prokka.annotation.sh <PROKKAGFF> [PFAM,TIGR,COG,EC]"
    exit 0
fi

if [ "$type" == "PFAM" ] ; then
    egrep -w "Dbxref=pfam[0-9]{5}" $infile | cut -f9| sed "s/^ID=\(PROKKA_[MOD_]*[0-9]\+\).\+Dbxref=pfam\([0-9]\+\).\+/\1\tPF\2/g"
elif [ "$type" == "TIGR" ] ; then
    egrep -w "Dbxref=TIGR[0-9]{5}" $infile | cut -f9| sed "s/^ID=\(PROKKA_[MOD_]*[0-9]\+\).\+Dbxref=TIGR\([0-9]\+\).\+/\1\tTIGR\2/g"
elif [ "$type" == "EC" ] ; then
    grep "eC_number=" $infile | cut -f9 | cut -f1,2 -d ';'| sed 's/ID=//g'| sed 's/;eC_number=/\t/g'
else
    egrep "Dbxref=COG[0-9]{4}" $infile | cut -f9 | sed 's/^ID=\(PROKKA_[MOD_]*[0-9]\+\).\+COG\([0-9]\+\),.\+/\1\tCOG\2/g'
fi
