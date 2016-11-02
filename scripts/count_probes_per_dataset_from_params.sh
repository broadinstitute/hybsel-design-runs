#!/bin/bash

# Input:
#  - $1: a parameters file output by find_optimal_params.py
#  - $2: a directory containing folders that correspond to each
#        dataset, each of which contains a fasta file for various
#        parameters that give probe sequences (e.g.,
#        $2/lassa/[parameters].fasta gives probe sequences for the
#        'lassa' dataset with [parameters])
#
# Output to stdout:
#  each dataset with the number of probes from that dataset used in
#  the design given by the parameters in $1

while read line; do
    dataset=$(echo "$line" | awk -F'\t' '{print $1}')
    params=$(echo "$line" | awk -F'\t' '{print $2}')
    params_sep=$(echo "$params" | sed 's/(//' | sed 's/)//' | sed 's/, / /')
    mismatches=$(echo "$params_sep" | awk '{print $1}')
    coverextension=$(echo "$params_sep" | awk '{print $2}')
    fasta="$2/$dataset/mismatches_${mismatches}-coverextension_${coverextension}.n_expanded.fasta"
    num_probes=$(cat $fasta | grep -v '>' | wc -l)
    echo -e "$dataset\t$num_probes"
done < $1
