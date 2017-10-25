#!/bin/bash

# Convert a result set from the Influenza Virus Resource to a format
# identical to NCBI's viral accession list.

echo -n "" > influenza.acc-list.txt
while read line; do
    acc=$(echo "$line" | awk -F'\t' '{print $1}')
    segment=$(echo "$line" | awk -F'\t' '{split($4, s, " "); print s[1]}')
    name=$(echo "$line" | awk -F'\t' '{print $9}')
    flu_type=$(echo "$name" | awk -F' ' '{print $2}')

    refseq=$(awk -F'\t' -v t="Influenza $flu_type virus" -v s="segment $segment" '$1==t && $2==s {print $3}' influenza.refseq.txt)

    echo -e "$refseq\t$acc\thuman\tOrthomyxoviridae,Influenzavirus $flu_type,Influenza $flu_type virus\t$name\tsegment $segment" >> influenza.acc-list.txt
done < <(grep -v '^accession' influenza.virus-resource-download.txt)
