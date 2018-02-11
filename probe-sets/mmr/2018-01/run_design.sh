#!/bin/bash

# Args:
#   1: subcommand
#
#   If subcommand is "make-design-commands":
#      2: comma-separated list of mismatches to search
#      3: comma-separated list of cover extensions to search
#   If subcommand is "parallel":
#      2: number of jobs to run in parallel

mkdir -p designs

DESIGN_COMMANDS="commands.design.txt"
NUM_PROBES_TSV="designs/num-probes.tsv"

if [[ $1 == "make-design-commands" ]]; then
    # Make commands to call CATCH over a grid of parameter values

    IFS=',' read -ra MISMATCHES_TO_SEARCH <<< "$2"
    IFS=',' read -ra EXTENSIONS_TO_SEARCH <<< "$3"

    echo -n "" > $DESIGN_COMMANDS
    mkdir -p designs/probes
    mkdir -p designs/log

    for dataset in measles mumps rubella; do
        for m in "${MISMATCHES_TO_SEARCH[@]}"; do
            for e in "${EXTENSIONS_TO_SEARCH[@]}"; do
                outfnprefix="${dataset}.m${m}-e${e}"

                if ! [[ -f "$outfn" && -s "$outfn" ]]; then
                    cmd="design.py $dataset additional-input/${dataset}.fasta -pl 100 -ps 50 -l 100 -m $m -e $e --add-adapters --add-reverse-complements --expand-n -o designs/probes/${outfnprefix}.fasta --max-num-processes 8 --verbose &> designs/log/${outfnprefix}.err"
                    echo "$cmd" >> $DESIGN_COMMANDS
                fi
            done
        done
    done
elif [[ $1 == "parallel-design" ]]; then
    # Run commands to calculate number of probes

    NJOBS="$2"

    source activate /ebs/hybsel-design-runs/tools/envs/catch-prod

    parallel --jobs $NJOBS --no-notice --progress < $DESIGN_COMMANDS

    source deactivate
elif [[ $1 == "aggregate-probe-counts" ]]; then
    # Combine the number of probes for each point on the grid of
    # parameter values into one file; this file can be used to
    # pool probes across datasets

    echo -e "dataset\tmismatches\tcover_extension\tnum_probes" > $NUM_PROBES_TSV

    for f in $(ls -1 designs/probes); do
        numprobes=$(cat designs/probes/$f | grep '>' | wc -l)
        dataset=$(echo "$f" | awk -F'.' '{print $1}')

        params=$(echo "$f" | awk -F'.' '{print $2}')
        m=$(echo "$params" | awk -F'-' '{print $1}' | sed 's/m//')
        e=$(echo "$params" | awk -F'-' '{print $2}' | sed 's/e//')

        echo -e "$dataset\t$m\t$e\t$numprobes" >> $NUM_PROBES_TSV
    done
else
    echo "FATAL: Unknown subcommand"
    exit 1
fi
