#!/bin/bash

COMMANDS="commands.txt"

if [[ $1 == "make-commands" ]]; then
    echo -n "" > $COMMANDS

    # Generate a space-separated list of genomes to blacklist; put
    # 'custom:' in front of each so that make_probes.py will read
    # it as a custom fasta input
    blacklisted_genomes=$(ls -1 blacklist-genomes-subset/*.fa | awk '{print "custom:"$0}' | tr '\n' ' ' | sed 's/ $//')

    # Use '--small_seq_min 60' because a small number of sequences (1475)
    # have a length that is less than the probe length, and the length of the
    # shortest sequence is 60.
    for m in 2 4 6; do
        for e in 25 50 100; do
            mkdir -p probes/m${m}-e${e}-iem30-expandn
            for x in input/*; do
                f=$(echo "$x" | awk -F'/' '{print $2}')
                cmd="make_probes.py -pl 75 -ps 25 -l 75 -m $m -e $e --island_of_exact_match 30 -mt 8 --island_of_exact_match_tolerant 0 --expand_n --small_seq_min 60 --skip_reverse_complements --skip_adapters -d custom:${x} --blacklist_genomes $blacklisted_genomes -o probes/m${m}-e${e}-iem30-expandn/${f}.fasta --use_native_dict_when_finding_tolerant_coverage --max_num_processes 2 --verbose > log/m${m}-e${e}-iem30-expandn.${f}.out 2> log/m${m}-e${e}-iem30-expandn.${f}.err"
                echo "$cmd" >> $COMMANDS
            done
        done
    done
elif [[ $1 == "parallel" ]]; then
    source activate /ebs/hybsel-design-runs/tools/envs/hybseldesign-prod

    parallel --jobs 20 --no-notice --progress < $COMMANDS

    source deactivate
else
    echo "FATAL: Unknown subcommand"
    exit 1
fi
