#!/bin/bash

#for d in dengue_1 dengue_2 dengue_3 dengue_4 zika chikungunya; do
#for d in dengue_1 dengue_2 dengue_3 dengue_4; do
for d in west_nile; do
    for m in 0 1 2 3 4 5; do
        for e in 0 10 20 30 40 50; do
            if [[ "$d" =~ ^dengue* ]]; then
                mem="32g"
            elif [[ "$d" = "west_nile" ]]; then
                mem="40g"
            else
                mem="8g"
            fi
            echo "Submitting $d (m=$m, e=$e, mem=$mem)"
            qbsub --profile -l h_vmem=${mem} -q long -N ${d}-${m}-${e} -o ${d}/mismatches_${m}-coverextension_${e}.out python ~/viral/hybsel_design/bin/make_probes.py -pl 100 -ps 50 -c 1.0 -d $d -m $m -e $e -o ${d}/mismatches_${m}-coverextension_${e}.fasta --print_analysis --verbose
        done
    done
done
