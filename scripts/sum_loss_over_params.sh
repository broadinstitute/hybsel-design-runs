#!/bin/bash

# Given a parameters file output by find_optimal_params.py, compute
# a loss of the parameters by computing a loss for each dataset
# (virus) and then summing over these losses.

# The loss for a single dataset is defined by the 'mismatches' and
# 'cover_extension' parameter values assigned to that dataset:
#  LOSS = mismatches^2 + (cover_extension / 10.0)^2

awk -F'\t' '{print $2}' $1 | sed 's/(//' | sed 's/)//' | sed 's/, /\t/' | awk '{print $1*$1 + ($2/10.0)*($2/10.0)}' | awk '{sum+=$0} END {print sum}'
