#!/bin/bash

# Args:
#   1: subcommand
#   2: dataset to generate data for
#
#   If subcommand is "make-commands":
#      3: sample sizes (i.e., number of genomes) to use
#         as input; 3 comma-separated values: start,end,stride
#         where start is the smallest sample size, end is
#         the largest, and stride is the space between them
#      4: approaches to try (comma separated) - 'sc', 'nt',
#         'nr', and 'ds'
#      5: mismatches to try for approaches where it is
#         applicable (comma separated)
#      6: cover extensions to try for set cover approach
#         (comma separated)

DATASET=$2
mkdir -p $DATASET

# Draw samples of the following sizes: between start and
# end, spaced by the stride; the sample size is the number
# of genomes randomly selected to be used as an input (i.e.,
# in the plot generated downstream, it is the number of
# genomes as plotted on the x-axis)
SAMPLE_SIZE_START=$(echo "$3" | awk -F, '{print $1}')
SAMPLE_SIZE_END=$(echo "$3" | awk -F, '{print $2}')
SAMPLE_SIZE_STRIDE=$(echo "$3" | awk -F, '{print $3}')

# Draw the following number of samples for each value of
# the sample size (used to construct a confidence interval
# for a particular combination of approach and sample
# size)
NUM_SAMPLES=5

IFS=',' read -ra APPROACHES_TO_TRY <<< "$4"
IFS=',' read -ra MISMATCHES_TO_TRY <<< "$5"
IFS=',' read -ra EXTENSIONS_TO_TRY <<< "$6"

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

COMMANDS="$DATASET/commands.txt"
COMBINED_NUM_PROBES="$DATASET/num-probes.tsv"

if [[ $1 == "make-commands" ]]; then
    # Make commands to call hybsel_design to calculate number of probes

    echo -n "" > $COMMANDS
    mkdir -p $DATASET/num-probes
    mkdir -p $DATASET/log

    for n in $(seq $SAMPLE_SIZE_START $SAMPLE_SIZE_STRIDE $SAMPLE_SIZE_END); do
        for k in $(seq 1 $NUM_SAMPLES); do

            # Design probes with set cover approach
            if try_approach "sc"; then
                for m in "${MISMATCHES_TO_TRY[@]}"; do
                    for e in "${EXTENSIONS_TO_TRY[@]}"; do
                        outfn="$DATASET/num-probes/sc.m${m}-e${e}.n${n}.k${k}"
                        if ! [[ -f "$outfn" && -s "$outfn" ]]; then
                            cmd="make_probes.py -pl 75 -ps 25 -l 75 -m $m -e $e --skip_reverse_complements --skip_adapters -d $DATASET --limit_target_genomes_randomly_with_replacement $n --max_num_processes 8 --verbose > $outfn 2> $DATASET/log/sc.m${m}-e${e}.n${n}.k${k}.err"
                            echo "$cmd" >> $COMMANDS
                        fi
                    done
                done
            fi

            # Design probes with naive tiling approach
            if try_approach "nt"; then
                outfn="$DATASET/num-probes/nt.n${n}.k${k}"
                if ! [[ -f "$outfn" && -s "$outfn" ]]; then
                    cmd="make_probes_naively.py -pl 75 -ps 25 --skip_reverse_complements -d $DATASET --limit_target_genomes_randomly_with_replacement $n --verbose > $outfn 2> $DATASET/log/nt.n${n}.k${k}.err"
                    echo "$cmd" >> $COMMANDS
                fi
            fi

            # Design probes with naive redundant filter
            if try_approach "nr"; then
                for m in "${MISMATCHES_TO_TRY[@]}"; do
                    outfn="$DATASET/num-probes/nr.m${m}.n${n}.k${k}"
                    if ! [[ -f "$outfn" && -s "$outfn" ]]; then
                        cmd="make_probes_naively.py -pl 75 -ps 25 --naive_redundant_filter $m 65 --skip_reverse_complements -d $DATASET --limit_target_genomes_randomly_with_replacement $n --verbose > $outfn 2> $DATASET/log/nr.m${m}.n${n}.k${k}.err"
                        echo "$cmd" >> $COMMANDS
                    fi
                done
            fi

            # Design probes with dominating set filter
            if try_approach "ds"; then
                for m in "${MISMATCHES_TO_TRY[@]}"; do
                    outfn="$DATASET/num-probes/ds.m${m}.n${n}.k${k}"
                    if ! [[ -f "$outfn" && -s "$outfn" ]]; then
                        cmd="make_probes_naively.py -pl 75 -ps 25 --dominating_set_filter $m 65 --skip_reverse_complements -d $DATASET --limit_target_genomes_randomly_with_replacement $n --verbose > $outfn 2> $DATASET/log/ds.m${m}.n${n}.k${k}.err"
                        echo "$cmd" >> $COMMANDS
                    fi
                done
            fi

        done
    done
elif [[ $1 == "parallel" ]]; then
    # Run commands to calculate number of probes

    source activate /ebs/hybsel-design-runs/tools/envs/hybseldesign-prod

    parallel --jobs 36 --no-notice --progress < $COMMANDS

    source deactivate
elif [[ $1 == "merge-counts" ]]; then
    # Combine all output files, each of which has the number of probes from
    # a run of hybsel_design, into one file; this single file can be used
    # for producing a plot

    echo -n "" > $COMBINED_NUM_PROBES

    for f in $(ls -1 $DATASET/num-probes); do
        num=$(cat $DATASET/num-probes/$f | tr -d '\n')
        approach=$(echo "$f" | awk -F'.' '{print $1}')

        if [[ "$approach" == "sc" ]]; then
            params=$(echo "$f" | awk -F'.' '{print $2}')
            m=$(echo "$params" | awk -F'-' '{print $1}' | sed 's/m//')
            e=$(echo "$params" | awk -F'-' '{print $2}' | sed 's/e//')
            n=$(echo "$f" | awk -F'.' '{print $3}' | sed 's/n//')

            echo -e "Set cover, m=$m, e=$e\t$n\t$num" >> $COMBINED_NUM_PROBES
        elif [[ "$approach" == "nt" ]]; then
            n=$(echo "$f" | awk -F'.' '{print $2}' | sed 's/n//')

            echo -e "Tiling\t$n\t$num" >> $COMBINED_NUM_PROBES
        elif [[ "$approach" == "nr" ]]; then
            m=$(echo "$f" | awk -F'.' '{print $2}' | sed 's/m//')
            n=$(echo "$f" | awk -F'.' '{print $3}' | sed 's/n//')

            echo -e "Naive redundant, m=$m\t$n\t$num" >> $COMBINED_NUM_PROBES
        elif [[ "$approach" == "ds" ]]; then
            m=$(echo "$f" | awk -F'.' '{print $2}' | sed 's/m//')
            n=$(echo "$f" | awk -F'.' '{print $3}' | sed 's/n//')

            echo -e "Dominating set, m=$m\t$n\t$num" >> $COMBINED_NUM_PROBES
        else
            echo "FATAL: Unknown approach for file $f in $DATASET/num-probes/"
            exit 1
        fi
    done
else
    echo "FATAL: Unknown subcommand"
    exit 1
fi
