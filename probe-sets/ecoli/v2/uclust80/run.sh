#!/bin/bash

COMMANDS="commands.txt"

if [[ $1 == "make-commands" ]]; then
    echo -n "" > $COMMANDS

    # Use '--small_seq_min 60' because a small number of sequences (1475)
    # have a length that is less than the probe length, and the length of the
    # shortest sequence is 60.
    for m in 2; do
        for e in 25; do
            mkdir -p probes/m${m}-e${e}-iem30
            for x in input/*; do
                f=$(echo "$x" | awk -F'/' '{print $2}')
                cmd="make_probes.py -pl 75 -ps 25 -l 75 -m $m -e $e --island_of_exact_match 30 --small_seq_min 60 --skip_reverse_complements --skip_adapters -d custom:${x} -o probes/m${m}-e${e}-iem30/${f}.fasta --max_num_processes 2 --verbose > log/m${m}-e${e}-iem30.${f}.out 2> log/m${m}-e${e}-iem30.${f}.err"
                echo "$cmd" >> $COMMANDS
            done
        done
    done
elif [[ $1 == "parallel" ]]; then
    source activate /ebs/hybsel-design-runs/tools/envs/hybseldesign-prod

    parallel --no-notice --progress < $COMMANDS

    source deactivate
else
    echo "FATAL: Unknown subcommand"
    exit 1
fi
