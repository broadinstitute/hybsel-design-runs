#!/bin/bash

# Design for a 12,472 array (synthesized by CustomArray)
# Use '--round-params 1 5' because we searched over a grid with all mismatch
#  values and with cover extension values spaced by 5
# Use '--loss-coeffs 1 0.0025' so that a cover extension of 20 has the same
#  loss as 1 mismatch (i.e., 1 * 1^2 = 0.0025 * 20^2) (similarly, so that
#  a cover extension of 100 has the same loss as 5 mismatches -- i.e.,
#  1 * 5^2 = 0.0025 * 100^2); this choice seems reasonable

mkdir -p pool

pool.py designs/num-probes.tsv 12472 pool/param-choices.tsv --round-params 1 5 --loss-coeffs 1 0.0025 --verbose &> pool/pool.out
