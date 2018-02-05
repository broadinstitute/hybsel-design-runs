#!/bin/bash

# Args:
#   1: subcommand
#   2: dataset to generate data for
#
#   If subcommand is "make-commands":
#      3: sample size (i.e., number of genomes) to use as input
#         for each experiment
#      4: mismatches to try (comma separated)
#      5: cover extensions to try (comma separated)
#   If subcommand is "parallel":
#      3: number of jobs to run in parallel

DATASET=$2
mkdir -p $DATASET

COMMANDS="$DATASET/commands.txt"
COMBINED_NUM_PROBES="$DATASET/num-probes.tsv"

if [[ $1 == "make-commands" ]]; then
    # Make commands to call CATCH to calculate number of probes

    SAMPLE_SIZE="$3"

    # Draw the following number of samples for each tested point
    # (i.e., choice of mismatches and cover extension); this is
    # used to construct a confidence interval at each point
    NUM_SAMPLES=5

    IFS=',' read -ra MISMATCHES_TO_TRY <<< "$4"
    IFS=',' read -ra EXTENSIONS_TO_TRY <<< "$5"

    # Check if array contains an element
    # From https://stackoverflow.com/a/8574392
    contains_element() {
      local e match="$1"
      shift
      for e; do [[ "$e" == "$match" ]] && return 0; done
      return 1
    }
    try_approach() {
        return $(contains_element $1 "${APPROACHES_TO_TRY[@]}")
    }

    echo -n "" > $COMMANDS
    mkdir -p $DATASET/num-probes
    mkdir -p $DATASET/log

    for k in $(seq 1 $NUM_SAMPLES); do
        for m in "${MISMATCHES_TO_TRY[@]}"; do
            for e in "${EXTENSIONS_TO_TRY[@]}"; do
                outfn="$DATASET/num-probes/sc.m${m}-e${e}.k${k}"

                if ! [[ -f "$outfn" && -s "$outfn" ]]; then
                    cmd="design.py $DATASET -pl 75 -ps 25 -l 75 -m $m -e $e --limit-target-genomes-randomly-with-replacement $SAMPLE_SIZE --max-num-processes 8 --verbose > $outfn 2> $DATASET/log/sc.m${m}-e${e}.k${k}.err"
                    echo "$cmd" >> $COMMANDS
                fi
            done
        done
    done
elif [[ $1 == "parallel" ]]; then
    # Run commands to calculate number of probes

    NJOBS="$3"

    source activate /ebs/hybsel-design-runs/tools/envs/catch-prod

    parallel --jobs $NJOBS --no-notice --progress < $COMMANDS

    source deactivate
elif [[ $1 == "merge-counts" ]]; then
    # Combine all output files, each of which has the number of probes from
    # a run of CATCH, into one file; this single file can be used
    # for producing a plot

    echo -e "approach\tmismatches\textension\tnumprobes" > $COMBINED_NUM_PROBES

    for f in $(ls -1 $DATASET/num-probes); do
        num=$(cat $DATASET/num-probes/$f | tr -d '\n')
        approach=$(echo "$f" | awk -F'.' '{print $1}')

        params=$(echo "$f" | awk -F'.' '{print $2}')
        m=$(echo "$params" | awk -F'-' '{print $1}' | sed 's/m//')
        e=$(echo "$params" | awk -F'-' '{print $2}' | sed 's/e//')

        echo -e "$approach\t$m\t$e\t$num" >> $COMBINED_NUM_PROBES
    done
else
    echo "FATAL: Unknown subcommand"
    exit 1
fi
