#!/usr/bin/env bash
if [ "$infile" == "" ] ; then
    echo "Usage: classification_to_krona.sh <input classification csv> <output prefix>"
    exit 0
fi
cut -f 3- -d ',' $1 | sort | uniq -c | sort -n | sed 's/,/	/g' | sed 's/\([1-9]\) /\1	/' > $2.krona.tab
ktImportText -o $2.krona.html $2.krona.tab
