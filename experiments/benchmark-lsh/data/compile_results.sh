#!/bin/bash

OUTFN="stats.tsv"
echo -e "data\tsample.size\tlsh\tnum.probes\ttime\tmem" > $OUTFN

DATASETS=("hiv1_without_ltr" "hepatitis_c" "zika" "ebola_zaire_with_2014" "influenza_a_segment4")
MISMATCHES_TO_TRY=("2" "5")
EXTENSIONS_TO_TRY=("25")

for dataset in "${DATASETS[@]}"; do
    if [ "$dataset" == "hiv1_without_ltr" ]; then
        SAMPLE_SIZES=(100 1000)
    elif [ "$dataset" == "hepatitis_c" ]; then
        SAMPLE_SIZES=(100 1000)
    elif [ "$dataset" == "zika" ]; then
        SAMPLE_SIZES=(100 350)
    elif [ "$dataset" == "ebola_zaire_with_2014" ]; then
        SAMPLE_SIZES=(100 800)
    elif [ "$dataset" == "influenza_a_segment4" ]; then
        SAMPLE_SIZES=(100 1000 5000 10000)
    fi

    for m in "${MISMATCHES_TO_TRY[@]}"; do
        for e in "${EXTENSIONS_TO_TRY[@]}"; do
            for sample_size in "${SAMPLE_SIZES[@]}"; do
                # No LSH
                outfn_prefix="$dataset/n${sample_size}.m${m}.e${e}.no-lsh"
                num_probes=$(grep '>' ${outfn_prefix}.fasta | wc -l)
                time=$(tail -n 4 ${outfn_prefix}.err | grep 'real' | awk '{print $2}')
                mem=$(tail -n 4 ${outfn_prefix}.err | grep 'mrss' | awk '{print $2}')
                echo -e "$dataset,m=$m,e=$e\t$sample_size\tnolsh\t$num_probes\t$time\t$mem" >> $OUTFN

                # LSH w/ Hamming
                outfn_prefix="$dataset/n${sample_size}.m${m}.e${e}.lsh-hamming"
                num_probes=$(grep '>' ${outfn_prefix}.fasta | wc -l)
                time=$(tail -n 4 ${outfn_prefix}.err | grep 'real' | awk '{print $2}')
                mem=$(tail -n 4 ${outfn_prefix}.err | grep 'mrss' | awk '{print $2}')
                echo -e "$dataset,m=$m,e=$e\t$sample_size\thamming\t$num_probes\t$time\t$mem" >> $OUTFN

                # LSH w/ MinHash
                outfn_prefix="$dataset/n${sample_size}.m${m}.e${e}.lsh-minhash"
                num_probes=$(grep '>' ${outfn_prefix}.fasta | wc -l)
                time=$(tail -n 4 ${outfn_prefix}.err | grep 'real' | awk '{print $2}')
                mem=$(tail -n 4 ${outfn_prefix}.err | grep 'mrss' | awk '{print $2}')
                echo -e "$dataset,m=$m,e=$e\t$sample_size\tminhash\t$num_probes\t$time\t$mem" >> $OUTFN
            done
        done
    done
done
