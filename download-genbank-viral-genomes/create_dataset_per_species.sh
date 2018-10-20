#!/bin/bash

# Create a datasets file, consisting of one dataset per
# species, from a list of lineages.
#
# Args:
#   1: file giving list of lineages (one per row, with
#      family/genus/species tab-separated)
#
# This prints to stdout two columns: col 1 gives a dataset
# name and col 2 gives the species corresponding to that dataset.

dataset_name() {
    # Create a dataset name from a species name
    # Convert to lowercase, make spaces be underscores, and if it
    # ends in ' virus' then remove that
    dataset=$(echo "$1" | awk '{ print tolower($0) }' | sed 's/ /\_/g' | sed -E 's/\_virus$//')
    echo "$dataset"
}

# Write a file giving dataset (column 1) and species (column 2)
tmpf=$(mktemp)
while read -r lineage; do
    species=$(echo "$lineage" | awk -F'\t' '{print $3}')
    dataset=$(dataset_name "$species")
    echo -e "$dataset\t$species" >> $tmpf
done < "$1"

# Verify that the number of unique dataset names equals the
# number of lineages
num_lineages=$(cat "$1" | sort | uniq | wc -l)
num_datasets=$(cat "$tmpf" | sort | uniq | wc -l)
if [ "$num_lineages" -ne "$num_datasets" ]; then
    echo "FATAL: Number of unique dataset names ($num_datasets) does not match number of lineages ($num_lineages)"
    exit 1
fi

# Print out, and then remove, the datasets file
cat "$tmpf"

rm "$tmpf"
