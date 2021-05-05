#!/bin/bash

for ps in V-All.201606 V-All.350k.201810; do
    for m in 5 10; do
        for e in 25 50; do
            Rscript plot.R out/${ps}.m${m}.e${e}.tsv.gz plots/${ps}.m${m}.e${e}.pdf
        done
    done
done
