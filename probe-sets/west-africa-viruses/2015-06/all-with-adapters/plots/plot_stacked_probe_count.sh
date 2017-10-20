#!/bin/bash

# Leave out final argument to label all datasets
Rscript ../../../../../scripts/plot_stacked_probe_count_by_dataset.R ../params.90000.probe_count_per_dataset.txt params.90000.probe_count_per_dataset.pdf 0,90000,30000
