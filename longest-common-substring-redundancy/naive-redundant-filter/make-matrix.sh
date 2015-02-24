#!/bin/bash
echo "" > matrix.vals
for l in $(seq 70 100); do
  for m in $(seq 0 12); do
    bsub -o make-matrix.out -q week -R "rusage[mem=1]" "python ~/viral/hybsel_design/bin/make_ebola_zaire_probes_by_lcs.py -l $l -m $m >> matrix.vals"
  done
done
