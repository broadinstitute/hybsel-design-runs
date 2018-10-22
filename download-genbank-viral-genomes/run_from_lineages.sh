#!/bin/bash

# Download sequences, and create other files, all from a
# list of lineages to include and (modified) version of
# the viral accession list.
# Note that Influenza should be handled separately (see
# influenza/ dir).
# Note also that it may be useful to clear the all/
# and accession_num/ directories before running this.

LINEAGES="lineages.txt"
ACC_LIST="taxid10239.nbr.20181019.modified.txt"
HUMAN_HOST_LINEAGES_TO_EXPLICITLY_INCLUDE="human-host-lineages-to-add.txt"

DATASETS_TO_CREATE="datasets.txt"
SEQS_LIST_TO_CREATE="seqs-in-201810-pull.txt"
DATASET_LINEAGES_TO_CREATE="seqs-in-201810-pull.lineages-per-dataset.txt"


# Create a datasets file
./create_dataset_per_species.sh $LINEAGES > $DATASETS_TO_CREATE

# Run the download_dataset_fastas.py script to make sure everything
# (sequences and datasets) matches up, without doing anything; it will
# raise errors if not
python download_dataset_fastas.py -dl $DATASETS_TO_CREATE -gl $ACC_LIST --human-host-lineages-to-add $HUMAN_HOST_LINEAGES_TO_EXPLICITLY_INCLUDE --skip-download -o all

# Print a list of sequences to include
python download_dataset_fastas.py -dl $DATASETS_TO_CREATE -gl $ACC_LIST --human-host-lineages-to-add $HUMAN_HOST_LINEAGES_TO_EXPLICITLY_INCLUDE --skip-download -o all --print-sequences > $SEQS_LIST_TO_CREATE

# Create a list of lineages per dataset (since the datasets file
# was created with create_dataset_per_species.sh, this should only
# have one lineage per dataset)
./convert_seqs_pull_to_lineages.sh $SEQS_LIST_TO_CREATE > $DATASET_LINEAGES_TO_CREATE

# Write a list of accession numbers for each dataset
python download_dataset_fastas.py -dl $DATASETS_TO_CREATE -gl $ACC_LIST --human-host-lineages-to-add $HUMAN_HOST_LINEAGES_TO_EXPLICITLY_INCLUDE --skip-download -o all --write-accession-nums accession_nums

# Download the sequences and create .py files, putting them into all/
python download_dataset_fastas.py -dl $DATASETS_TO_CREATE -gl $ACC_LIST --human-host-lineages-to-add $HUMAN_HOST_LINEAGES_TO_EXPLICITLY_INCLUDE --consolidate-segmented-genomes-into-one-fasta --gzip-fastas -o all
