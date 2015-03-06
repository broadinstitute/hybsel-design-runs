#!/bin/bash
echo "" > matrix.limit-target-genome.vals
for ltg in $(seq 10 10 140); do
  bsub -o make-matrix.out -q forest -R "rusage[mem=1]" "python ~/viral/hybsel_design/bin/make_probes_by_redundancy_reduction.py -d ebola_zaire -l 100 -m 0 -f ds --limit_target_genomes $ltg >> matrix.limit-target-genomes.vals"
  bsub -o make-matrix.out -q forest -R "rusage[mem=1]" "python ~/viral/hybsel_design/bin/make_probes_by_redundancy_reduction.py -d ebola_zaire -l 85 -m 2 -f ds --limit_target_genomes $ltg >> matrix.limit-target-genomes.vals"
done
