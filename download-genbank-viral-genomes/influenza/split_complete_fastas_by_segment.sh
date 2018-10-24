#!/bin/bash

# The Influenza A/B/C FASTA files created by download_influenza_seqs.sh
# contain all segments together (i.e., one FASTA file per species; with
# all 8 segments in each FASTA). This splits each FASTA file into 8
# separate ones: one per segment.
#
# This assumes the datasets are named 'influenza_{a,b,c}' and the data
# is located in ../all/data/ and are gzip'd.

DATASET_PYTHON_TEMPLATE_UNSEGMENTED="../dataset_unsegmented.template.py"

for dataset in influenza_a influenza_b influenza_c; do
    dataset_fasta="../all/data/${dataset}.fasta"

    # Decompress the FASTA
    gzip -d ${dataset_fasta}.gz

    for segment in 1 2 3 4 5 6 7 8; do
        if [ "$dataset" == "influenza_c" ] && [ "$segment" == "8" ]; then
            # Influenza C virus has only 7 segments, so skip this
            continue
        fi

        # Make a FASTA file for this segment
        segment_fasta="../all/data/${dataset}_segment${segment}.fasta"

        # Use the split_by_regex.py script to extract all sequences for
        # this segment
        python ~/misc-scripts/fasta-processing/split_by_regex.py -i $dataset_fasta -r ".*\[segment $segment\].*" -o1 ${segment_fasta} -o2 /dev/null

        num_seqs=$(cat $segment_fasta | grep '>' | wc -l)
        species=$(awk -F'\t' -v d="$dataset" '$1==d {print $2}' influenza.datasets.txt)

        # Make a dataset .py file for this segment
        segment_dataset_py="../all/${dataset}_segment${segment}.py"
        cat $DATASET_PYTHON_TEMPLATE_UNSEGMENTED | sed "s/\[\[VIRUS_REGEX\]\]/$species [segment $segment]/g" | sed "s/\[\[NUM_GENOMES\]\]/$num_seqs/g" | sed "s/\[\[SUBSET_NOTE\]\]//g" | sed "s/\[\[DATASET_NAME\]\]/${dataset}_segment${segment}/g" | sed "s/\[\[GZIP\]\]/.gz/g" > $segment_dataset_py

        # Compress the FASTA
        gzip ${segment_fasta}
    done

    # Recompress the FASTA
    gzip ${dataset_fasta}
done

echo "NOTE: In order for the influenza_{a,b,c}.py dataset files to use the combined segment data (instead of a separate data file with all segments), remember to edit those .py files!"
