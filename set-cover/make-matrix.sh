#!/bin/bash
#echo "" > matrix.vals
#for c in $(seq 0.1 0.1 1.0); do
#  bsub -o make-matrix.out -q forest -R "rusage[mem=2]" "python ~/viral/hybsel_design/bin/make_probes_by_set_cover.py -d ebola_zaire -l 100 -m 0 -c $c >> matrix.vals"
#done
#for c in $(seq 0.1 0.1 1.0); do
#  bsub -o make-matrix.out -q forest -R "rusage[mem=2]" "python ~/viral/hybsel_design/bin/make_probes_by_set_cover.py -d ebola_zaire -l 85 -m 2 -c $c >> matrix.vals"
#done
for l in $(seq 80 5 100); do
  for m in $(seq 0 1 5); do
    bsub -o make-matrix.out -q forest -R "rusage[mem=2]" "python ~/viral/hybsel_design/bin/make_probes_by_set_cover.py -d ebola_zaire -l $l -m $m -c 1.0 >> matrix.vals"
  done
done
