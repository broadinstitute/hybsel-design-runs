#!/bin/bash

# Input:
#  - $1: a parameters file output by find_optimal_params.py
#  - $2: a directory containing folders that correspond to each
#        dataset, each of which contains a fasta file for various
#        parameters that give probe sequences (e.g.,
#        $2/lassa/[parameters].fasta gives probe sequences for the
#        'lassa' dataset with [parameters])
#  - $3: reads as bam file
#  - $4: tmp folder
#
# Output to stdout:
#  each dataset that has at least one probe with an alignment to a
#  read, along with the number of reads that align to probes belonging
#  to that dataset

randchars=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)
tmp="$4/tmp-$randchars"
mkdir $tmp

# Loop over all datasets: put the concatenation of all the probes into
# $tmp/all_probes.fasta and put into $tmp/probe_datasets.txt each probe
# id along with the dataset that probe belongs to
echo -n "" > $tmp/all_probes.fasta
echo -n "" > $tmp/probe_datasets.txt
while read line; do
    dataset=$(echo "$line" | awk -F'\t' '{print $1}')
    params=$(echo "$line" | awk -F'\t' '{print $2}')
    params_sep=$(echo "$params" | sed 's/(//' | sed 's/)//' | sed 's/, / /')
    mismatches=$(echo "$params_sep" | awk '{print $1}')
    coverextension=$(echo "$params_sep" | awk '{print $2}')
    fasta="$2/$dataset/mismatches_${mismatches}-coverextension_${coverextension}.fasta"
    cat $fasta >> $tmp/all_probes.fasta
    grep '>' $fasta | awk '{print $1"\t'$dataset'"}' | sed 's/>//' >> $tmp/probe_datasets.txt
done < $1

# Align reads to $tmp/all_probes.fasta
align_reads_to_fasta $3 $tmp/all_probes.fasta $tmp &> $tmp/align.out

# Extract probes that a read aligns to
grep probe $tmp/alignment.sam | grep -v '@' | awk '{print $3}' | sort | uniq -c > $tmp/probes_with_alignment.uniq.txt

# Using $tmp/probe_datasets.txt, lookup the dataset for each probe 
echo -n "" > $tmp/probes_with_alignment.with_dataset.txt
while read probe_with_count; do
    arr=($probe_with_count)
    count=${arr[0]}
    probe=${arr[1]}
    dataset=$(grep "$probe" $tmp/probe_datasets.txt | awk '{print $2}')
    echo -e "$probe\t$count\t$dataset" >> $tmp/probes_with_alignment.with_dataset.txt
done < $tmp/probes_with_alignment.uniq.txt

# For each unique dataset, sum up the counts of the number of reads aligning
# to probes in that dataset
cat $tmp/probes_with_alignment.with_dataset.txt | awk '{sum[$3]+=$2} END {for (dataset in sum) {print dataset"\t"sum[dataset]}}' > $tmp/dataset_counts.txt

cat $tmp/dataset_counts.txt | sort -k2nr

# Cleanup
rm -rf $tmp
