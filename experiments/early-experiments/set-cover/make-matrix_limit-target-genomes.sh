#!/bin/bash
echo "" > matrix.limit-target-genomes.vals
for ltg in $(seq 10 10 140); do
  bsub -o make-matrix.out -q forest -R "rusage[mem=2]" "python ~/viral/hybsel_design/bin/make_probes_by_set_cover.py -d ebola_zaire -l 100 -m 0 -c 1.0 --limit_target_genomes $ltg >> matrix.limit-target-genomes.vals"
  bsub -o make-matrix.out -q forest -R "rusage[mem=2]" "python ~/viral/hybsel_design/bin/make_probes_by_set_cover.py -d ebola_zaire -l 85 -m 2 -c 1.0 --limit_target_genomes $ltg >> matrix.limit-target-genomes.vals"
done
