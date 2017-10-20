#!/bin/bash

# Label the following datasets
LABELS="hepatitis_c,hiv1_without_ltr,enterovirus_b,enterovirus_c,herpesvirus_5,influenza_a,gbv_c,lassa,dengue,zika,mumps"

Rscript ../../../../../scripts/plot_stacked_probe_count_by_dataset.R ../params_decided.probe_count_per_dataset.txt params_decided.probe_count_per_dataset.pdf 0,350000,50000 $LABELS
