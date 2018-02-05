#!/bin/bash

for dataset in ebola_zaire_with_2014 hepatitis_c hiv1_without_ltr zika; do
    Rscript make_plot.R data/${dataset}/num-probes.tsv plots/${dataset}.pdf $dataset
done
