#!/bin/bash
echo "" > matrix.vals
for s in $(seq 0 20); do
  for m in $(seq 0 12); do
    bsub -o make-matrix.out -q week -R "rusage[mem=1]" "python ~/viral/hybsel_design/bin/replicate_ebola_zaire_probes.py --shift $s --mismatch_thres $m >> matrix.vals"
  done
done
