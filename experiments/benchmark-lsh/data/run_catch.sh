#!/bin/bash


DATASETS=("hiv1_without_ltr" "hepatitis_c" "zika" "ebola_zaire_with_2014" "influenza_a_segment4")
MISMATCHES_TO_TRY=("2" "5")
EXTENSIONS_TO_TRY=("25")

source activate catch

for dataset in "${DATASETS[@]}"; do
    mkdir -p $dataset

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

    for sample_size in "${SAMPLE_SIZES[@]}"; do
        for m in "${MISMATCHES_TO_TRY[@]}"; do
            for e in "${EXTENSIONS_TO_TRY[@]}"; do
                echo "dataset=$dataset, n=$sample_size, m=$m, e=$e, no LSH"
                outfn_prefix="$dataset/n${sample_size}.m${m}.e${e}.no-lsh"
                /usr/bin/time -f 'real %e\nuser %U\nsys %S\nmrss %M' design.py $dataset -pl 75 -ps 25 -l 75 -m $m -e $e --limit-target-genomes-randomly-with-replacement $sample_size --max-num-processes 8 --verbose -o $outfn_prefix.fasta --write-analysis-to-tsv $outfn_prefix.tsv > $outfn_prefix.out 2> $outfn_prefix.err

                echo "dataset=$dataset, n=$sample_size, m=$m, e=$e, LSH w/ Hamming"
                outfn_prefix="$dataset/n${sample_size}.m${m}.e${e}.lsh-hamming"
                /usr/bin/time -f 'real %e\nuser %U\nsys %S\nmrss %M' design.py $dataset -pl 75 -ps 25 -l 75 -m $m -e $e --filter-with-lsh-hamming $m --limit-target-genomes-randomly-with-replacement $sample_size --max-num-processes 8 --verbose -o $outfn_prefix.fasta --write-analysis-to-tsv $outfn_prefix.tsv > $outfn_prefix.out 2> $outfn_prefix.err

                echo "dataset=$dataset, n=$sample_size, m=$m, e=$e, LSH w/ MinHash"
                outfn_prefix="$dataset/n${sample_size}.m${m}.e${e}.lsh-minhash"
                /usr/bin/time -f 'real %e\nuser %U\nsys %S\nmrss %M' design.py $dataset -pl 75 -ps 25 -l 75 -m $m -e $e --filter-with-lsh-minhash 0.5 --limit-target-genomes-randomly-with-replacement $sample_size --max-num-processes 8 --verbose -o $outfn_prefix.fasta --write-analysis-to-tsv $outfn_prefix.tsv > $outfn_prefix.out 2> $outfn_prefix.err
            done
        done
    done
done

source deactivate
