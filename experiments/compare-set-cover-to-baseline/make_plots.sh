#!/bin/bash

for dataset in ebola_zaire_with_2014 hepatitis_c hiv1_without_ltr; do
    for scale in "linear" "log"; do
        Rscript make_plot.R data/${dataset}/num-probes.tsv plots/${dataset}.${scale}.pdf $dataset $scale
    done
done
