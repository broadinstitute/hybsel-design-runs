#!/bin/bash

# Convert genomes spread across many fastas into one fasta.
#
# For segmented genomes, the script download_dataset_fastas.py groups the
# sequences (one per segment) into genomes and creates one fasta file per
# genome that has all of the segments in it. In general, this is OK, but
# for viruses with a lot of sequence data (e.g., influenza_a) it can create
# an unwieldy number of sequences. This combines all of those fasta files
# into one, using the sequence header to denote an ID for the genome.
#
# Args:
#  1: directory that contains fastas to combine
#  2: output fasta

echo -n "" > $2

for f in $(find $1 -type f -name '*.fasta'); do
    g=$(echo "${f##*/}" | sed 's/\.fasta//')
    while read line; do
        if [[ "$line" == ">"* ]]; then
            echo "$line [genome $g]" >> $2
        else
            echo $line >> $2
        fi
    done < $f
done

