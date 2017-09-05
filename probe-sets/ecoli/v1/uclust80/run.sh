#!/bin/bash

COMMANDS="commands.txt"

if [[ $1 == "make-commands" ]]; then
    echo -n "" > $COMMANDS

    for x in input/*; do
        f=$(echo "$x" | awk -F'/' '{print $2}')
        cmd="make_probes.py -pl 75 -ps 25 -m 2 -e 25 --skip_reverse_complements --skip_adapters -d custom:${x} -o probes/${f}.fasta --verbose > log/${f}.out 2> log/${f}.err"
        echo "$cmd" >> $COMMANDS
    done
elif [[ $1 == "parallel" ]]; then
    parallel --no-notice --progress < $COMMANDS
else
    echo "FATAL: Unknown subcommand"
    exit 1
fi
