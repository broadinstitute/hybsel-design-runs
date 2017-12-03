#!/bin/bash

for n in 100 500 1000; do
    # This assumes that min_n_string_length has been set to be larger than the
    # length of a probe, so that candidate probes are created by tiling with
    # a fixed stride (and not reset when there is a stretch of Ns).
    echo "From aligned input, n=$n:"
    python ~/hybsel-design/bin/make_probes_naively.py -pl 75 -ps 25 -d custom:hiv1_without_ltr.n${n}.aligned.gap-to-N.fasta --skip_reverse_complements

    echo "From unaligned input, n=$n:"
    python ~/hybsel-design/bin/make_probes_naively.py -pl 75 -ps 25 -d hiv1_without_ltr --skip_reverse_complements --limit_target_genomes_randomly_with_replacement $n

    echo ""
done
