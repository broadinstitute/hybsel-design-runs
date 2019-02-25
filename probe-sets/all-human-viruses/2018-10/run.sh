#!/bin/bash

# Args:
#   1: subcommand
#
#   If subcommand is "param-exploration-parallel-run":
#      2: number of jobs to run in parallel
#   If subcommand is "find-optimal-pooled-params":
#      2: number of probes to design
#      3: number of jobs to run in parallel
#   If subcommand is "make-probe-set-from-optimal-pooled-params":
#      2: number of probes to design
#
# Author: Hayden Metsky


DATASETSFN="datasets.txt"
COMMANDS="commands.txt"

# Populate the search space
DATASETS=()
while read dataset; do
    DATASETS+=("$dataset")
done < $DATASETSFN

MISMATCHES_TO_TRY=()
for m in $(seq 0 8); do
    MISMATCHES_TO_TRY+=("$m")
done

EXTENSIONS_TO_TRY=()
for e in $(seq 0 10 50); do
    EXTENSIONS_TO_TRY+=("$e")
done

ISLAND_OF_EXACT_MATCHES_TO_TRY=(25)


if [[ $1 == "param-exploration-make-commands" ]]; then
    # Make commands to call CATCH for designing probes
    # at different choices of parameters
    #
    # This sorts commands by the total size (number of nucleotides) of
    # the dataset; the reason is that the longer executions should be
    # started earlier, and it uses the size of the dataset as a proxy
    # for how long it will take to run
    #
    # Input: none
    # Output: $COMMANDS file listing commands for call CATCH, each for
    #         a different dataset and choice of parameters; the commands
    #         are sorted by the number of sequences in the dataset

    # First write a file of commands in which the first column gives
    # the size of the dataset and the second gives the command
    echo -n "" > /tmp/catch-cmds

    for dataset in "${DATASETS[@]}"; do
        mkdir -p dataset-runs/$dataset/probes
        mkdir -p dataset-runs/$dataset/coverage-analysis
        mkdir -p dataset-runs/$dataset/log

        dataset_size=$(zcat /ebs/hybsel-design-runs/tools/catch-prod/catch/datasets/data/${dataset}.fasta.gz | grep -v '>' | wc -c)

        for m in "${MISMATCHES_TO_TRY[@]}"; do
            # --filter-with-lsh-hamming should be commensurate with, but not
            # greater than, m. Use m
            filter_with_lsh_hamming=$m

            for e in "${EXTENSIONS_TO_TRY[@]}"; do
                for i in "${ISLAND_OF_EXACT_MATCHES_TO_TRY[@]}"; do
                    outprobesfn="dataset-runs/$dataset/probes/m${m}-e${e}-i${i}.fasta"
                    outanalysisfn="dataset-runs/$dataset/coverage-analysis/m${m}-e${e}-i${i}.tsv"
                    outlogfn="dataset-runs/$dataset/log/m${m}-e${e}-i${i}.err"

                    if ! [[ -f "${outprobesfn}.gz" && -s "${outprobesfn}.gz" ]]; then
                        # Create a command that designs probes and writes them to a file
                        cmd="design.py $dataset -pl 75 -ps 25 -l 75 -m $m -e $e --island-of-exact-match $i --filter-with-lsh-hamming $filter_with_lsh_hamming --filter-polya 20 4 --expand-n 2 -o $outprobesfn --write-analysis-to-tsv $outanalysisfn --small-seq-skip 74 --max-num-processes 8 --verbose &> $outlogfn; gzip $outprobesfn; gzip $outanalysisfn; gzip -f $outlogfn"
                        echo -e "$dataset_size\t$cmd" >> /tmp/catch-cmds
                    fi
                done
            done
        done
    done

    # Sort all commands by the dataset size (first column), and write
    # to $COMMANDS the commands in sorted order (from high to low)
    sort -t$'\t' -k1nr /tmp/catch-cmds | awk -F'\t' '{print $2}' > $COMMANDS
    rm /tmp/catch-cmds

elif [[ $1 == "param-exploration-parallel-run" ]]; then
    # Run commands that run CATCH for designing probes
    # at different choices of parameters
    # Input: $COMMANDS file listing commands from `param-exploration-make-commands`
    # Output: Result of calling those commands - i.e., probes written to
    #         dataset-runs/[dataset]/probes/[param choices].fasta.gz

    if [ -z "$2" ]; then
        echo "FATAL: Unknown number of commands to run in parallel"
        exit 1
    fi

    NJOBS="$2"

    source activate /ebs/hybsel-design-runs/tools/envs/catch-prod

    parallel --jobs $NJOBS --no-notice --progress < $COMMANDS

    source deactivate

elif [[ $1 == "param-exploration-combine-counts" ]]; then
    # Combine probe counts to construct a table of the number of
    # probes for each choice of parameters for each dataset
    # Input: Probes given in dataset-runs/[dataset]/probes/[param choice].fasta.gz
    # Output: Table combining all the probe counts, across datasets and
    #         choices of parameters - written to pooled/probe-counts.tsv
    mkdir -p pooled/
    outfn="pooled/probe-counts.tsv"

    # Print a header
    echo -e "dataset\tmismatches\tcover_extension\tisland_of_exact_match\tnum_probes" > $outfn

    # Print the number of probes for each combination of parameters
    for dataset in "${DATASETS[@]}"; do
        for m in "${MISMATCHES_TO_TRY[@]}"; do
            for e in "${EXTENSIONS_TO_TRY[@]}"; do
                for i in "${ISLAND_OF_EXACT_MATCHES_TO_TRY[@]}"; do
                    outprobesfn="dataset-runs/$dataset/probes/m${m}-e${e}-i${i}.fasta.gz"
                    # Only print if the file (count) exists
                    if [[ -f "$outprobesfn" && -s "$outprobesfn" ]]; then
                        num_probes=$(zcat "$outprobesfn" | grep '>' | wc -l)
                        echo -e "$dataset\t$m\t$e\t$i\t$num_probes" >> $outfn
                    fi
                done
            done
        done
    done

elif [[ $1 == "find-optimal-pooled-params" ]]; then
    # Run CATCH (pool.py) to find the optimal choice of parameters for
    # each dataset, using the table generated by `param-exploration-combine-counts`
    # Input: pooled/probe-counts.tsv, the table generated by
    #        `params-exploration-combine-counts`; also, a number of probes
    #        for the design
    # Output: optimal choice of parameters for each dataset, such that the
    #         total number of probes is as desired - written to table at
    #         pooled/design-[num probes]/params.tsv

    if [ -z "$2" ]; then
        echo "FATAL: Unknown number of probes to design"
        exit 1
    fi
    num_probes="$2"

    if [ -z "$3" ]; then
        echo "FATAL: Unknown number of commands to run in parallel"
        exit 1
    fi
    NJOBS="$3"

    mkdir -p pooled/design-${num_probes}

    # Extract columns for mismatches and cover_extension (leave out island_of_exact_match)
    # Note: Here, we assume that island_of_exact_match is a constant (and determine
    # and store that constant in $iem)
    if [[ $(cat pooled/probe-counts.tsv | tail -n +2 | awk '{print $4}' | sort | uniq | wc -l) -gt "1" ]]; then
        # island_of_exact_match is not a constant
        echo "FATAL: island_of_exact_match must be a constant"
        exit 1
    fi
    iem=$(cat pooled/probe-counts.tsv | tail -n +2 | awk '{print $4}' | sort | uniq)
    awk '{print $1"\t"$2"\t"$3"\t"$5}' pooled/probe-counts.tsv > pooled/probe-counts.2params.tsv

    mkdir -p /tmp/opt-pooled-params
    NRUNS=64

    # Create commands to run pool.py 64 times (each has a random initial start)
    echo -n "" > $COMMANDS
    for i in $(seq 1 $NRUNS); do
        cmd="pool.py pooled/probe-counts.2params.tsv $num_probes /tmp/opt-pooled-params/${i}.params --round-params 1 10 --loss-coeffs 1 0.01 &> /tmp/opt-pooled-params/${i}.out"
        echo "$cmd" >> $COMMANDS
    done

    # Run the commands
    source activate /ebs/hybsel-design-runs/tools/envs/catch-prod
    parallel --jobs $NJOBS --no-notice --progress < $COMMANDS
    source deactivate

    # Find the overall minimum loss across the runs, and use the parameters for this run
    min_loss="999999999.0"
    min_i="0"
    for i in $(seq 1 $NRUNS); do
        loss=$(cat /tmp/opt-pooled-params/${i}.out | grep '^Loss:' | awk '{print $2}')
        if (( $(echo "$loss < $min_loss" | bc -l ) )); then
            min_loss="$loss"
            min_i="$i"
        fi
    done
    cp /tmp/opt-pooled-params/${min_i}.params pooled/design-${num_probes}/params.tsv
    cp /tmp/opt-pooled-params/${min_i}.out pooled/design-${num_probes}/pool.out
    rm -rf /tmp/opt-pooled-params

    # Add back in the column giving island_of_exact_match
    mv pooled/design-${num_probes}/params.tsv pooled/design-${num_probes}/params.tsv.tmp
    cat pooled/design-${num_probes}/params.tsv.tmp | head -n 1 | awk '{print $0"\tisland_of_exact_match"}' > pooled/design-${num_probes}/params.tsv
    cat pooled/design-${num_probes}/params.tsv.tmp | tail -n +2 | awk -v iem="$iem" '{print $1"\t"$2"\t"$3"\t"iem}' >> pooled/design-${num_probes}/params.tsv
    rm pooled/design-${num_probes}/params.tsv.tmp

elif [[ $1 == "make-probe-set-from-optimal-pooled-params" ]]; then
    # Construct a final probe set from the optimal parameter choices for
    # each dataset, as determined by `find-optimal-pooled-params`
    # Input: Table giving optimal choice of parameters for each dataset,
    #        stored in pooled/design-[num probes]/params.tsv; also, a
    #        number of probes for the design
    # Output: Probe set across all datasets - written to
    #         pooled/design-[num probes]/probes.fasta.gz
    #         Also, number of probes for each dataset - written to
    #         pooled/design-[num probes]/probe-counts.tsv
    #         Also, coverage analyses across all datasets - written to
    #         pooled/design-[num probes]/coverage-analysis.tsv.gz

    if [ -z "$2" ]; then
        echo "FATAL: Unknown number of probes to design"
        exit 1
    fi
    num_probes=$2

    if [[ ! -f "pooled/design-${num_probes}/params.tsv" ]]; then
        echo "FATAL: Unknown optimal pooled params"
        exit 1
    fi

    # Verify (using the header) that the order of columns (parameters) in
    # params.tsv is as expected
    paramsfn="pooled/design-${num_probes}/params.tsv"
    header=$(cat $paramsfn | head -n 1 | sed 's/\t/,/g')
    if [[ "$header" != "dataset,mismatches,cover_extension,island_of_exact_match" ]]; then
        echo "FATAL: Column order in pooled params is not the expected order"
        exit 1
    fi

    # Combine all the probes
    pooledprobesfn="pooled/design-${num_probes}/probes.fasta"
    echo -n "" > $pooledprobesfn
    while read line; do
        dataset=$(echo "$line" | awk '{print $1}')
        m=$(echo "$line" | awk '{print $2}')
        e=$(echo "$line" | awk '{print $3}')
        i=$(echo "$line" | awk '{print $4}')
        probesfn="dataset-runs/$dataset/probes/m${m}-e${e}-i${i}.fasta.gz"
        zcat $probesfn >> $pooledprobesfn
    done < <(cat $paramsfn | tail -n +2)
    gzip $pooledprobesfn

    # Count the number of probes for each dataset
    pooledprobecountsfn="pooled/design-${num_probes}/probe-counts.tsv"
    echo -e "dataset\tnum_probes" > $pooledprobecountsfn
    while read line; do
        dataset=$(echo "$line" | awk '{print $1}')
        m=$(echo "$line" | awk '{print $2}')
        e=$(echo "$line" | awk '{print $3}')
        i=$(echo "$line" | awk '{print $4}')
        probesfn="dataset-runs/$dataset/probes/m${m}-e${e}-i${i}.fasta.gz"
        num_probes_in_dataset=$(zcat "$probesfn" | grep '>' | wc -l)
        echo -e "$dataset\t$num_probes_in_dataset" >> $pooledprobecountsfn
    done < <(cat $paramsfn | tail -n +2)

    # Combine all the coverage analyses
    pooledanalysisfn="pooled/design-${num_probes}/coverage-analysis.tsv"
    # Get the header of the analysis tsv
    while read line; do
        dataset=$(echo "$line" | awk '{print $1}')
        m=$(echo "$line" | awk '{print $2}')
        e=$(echo "$line" | awk '{print $3}')
        i=$(echo "$line" | awk '{print $4}')
        # Get the header of the analysis tsv
        analysisfn="dataset-runs/$dataset/coverage-analysis/m${m}-e${e}-i${i}.tsv"
        zcat $analysisfn | head -n 1 > $pooledanalysisfn
    done < <(cat $paramsfn | tail -n +2 | head -n 1)
    # Concatenate all analysis tsvs, leaving out the header
    while read line; do
        dataset=$(echo "$line" | awk '{print $1}')
        m=$(echo "$line" | awk '{print $2}')
        e=$(echo "$line" | awk '{print $3}')
        i=$(echo "$line" | awk '{print $4}')
        analysisfn="dataset-runs/$dataset/coverage-analysis/m${m}-e${e}-i${i}.tsv"
        zcat $analysisfn | tail -n +2 >> $pooledanalysisfn
    done < <(cat $paramsfn | tail -n +2)
    gzip $pooledanalysisfn
else
    echo "FATAL: Unknown subcommand"
    exit 1
fi
