#!/bin/bash

# Given a parameters file output by find_optimal_params.py, count
# the frequency of each combination of parameters. This can be
# used to create a 2D bubble graph of the parameters.

awk -F'\t' '{print $2}' $1 | sed 's/(//' | sed 's/)//' | sed 's/, /\t/' | sort -k1,2 | uniq -c | awk '{print $2"\t"$3"\t"$1}'
