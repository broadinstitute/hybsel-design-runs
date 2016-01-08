#!/bin/bash

# Given a parameters file output by find_optimal_params.py, list
# each of the viruses sorted by their parameters, so that those
# with the highest parameters are listed at the end.

cat $1 | sed 's/(//' | sed 's/)//' | sed 's/, /\t/' | sort -k2,3 | column -t
