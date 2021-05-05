#!/bin/bash

COMMANDS="/tmp/commands.txt"
echo -n "" > $COMMANDS

for ps in V-All.201606 V-All.350k.201810; do
    for m in 5 10; do
        for e in 25 50; do
            echo -e "analyze_probe_coverage.py -d data/ncov.fasta.gz -f data/${ps}.fasta.gz -m $m -e $e  -l 75 --print-analysis --write-sliding-window-coverage out/${ps}.m${m}.e${e}.tsv --verbose &> out/${ps}.m${m}.e${e}.out; gzip out/${ps}.m${m}.e${e}.out; gzip out/${ps}.m${m}.e${e}.tsv" >> $COMMANDS
        done
    done
done

parallel --jobs 16 --no-notice --progress < $COMMANDS
