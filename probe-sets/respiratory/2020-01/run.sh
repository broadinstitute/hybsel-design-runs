#!/bin/bash

# Args:
#   1: subcommand
#   Subcommands:
#     - param-exploration-make-commands
#     - param-exploration-parallel-run
#     - param-exploration-combine-counts
#     - find-optimal-pooled-params
#     - make-probe-set-from-optimal-pooled-params
#
#   If subcommand is "param-exploration-parallel-run":
#      2: number of jobs to run in parallel
#   If subcommand is "find-optimal-pooled-params":
#      2: number of probes to design
#   If subcommand is "make-probe-set-from-optimal-pooled-params":
#      2: number of probes to design
#
# Author: Hayden Metsky


COMMANDS="commands.txt"

# Populate the search space
DATASETS=()
while read -r dataset; do
    DATASETS+=("$dataset")
done < input/datasets.txt

MISMATCHES_TO_TRY=(0 1 2 3)
for m in $(seq 4 2 8); do
    MISMATCHES_TO_TRY+=("$m")
done

EXTENSIONS_TO_TRY=()
for e in $(seq 0 25 75); do
    EXTENSIONS_TO_TRY+=("$e")
done

ISLAND_OF_EXACT_MATCHES_TO_TRY=()
for i in 25; do
    ISLAND_OF_EXACT_MATCHES_TO_TRY+=("$i")
done

NCBI_API_KEY="321e5004f2502f604ee2e8858a22b993d608"

if [[ -d "/ebs/tmpfs/tmp" ]]; then
    TMPDIR="/ebs/tmpfs/tmp"
else
    TMPDIR="/tmp"
fi
export TMPDIR

# Allow activating conda environments
source ~/anaconda3/etc/profile.d/conda.sh


dataset_load () {
    if [[ $1 == "sars" ]]; then
        # Load from file
        echo "input/SARS-all.20200218.fasta.gz"
    else
        taxid=$(cat input/taxid.tsv | awk -F'\t' -v d="$1" '$1==d {print $2}')
        echo "download:$taxid"
    fi
}


if [[ $1 == "param-exploration-make-commands" ]]; then
    # Make commands to call CATCH for determining number of probes
    # at different choices of parameters
    # Input: none
    # Output: $COMMANDS file listing commands for call CATCH, each for
    #         a different dataset and choice of parameters

    echo -n "" > $COMMANDS

    for dataset in "${DATASETS[@]}"; do
        mkdir -p dataset-runs/$dataset

        dataset_input=$(dataset_load $dataset)

        # For influenza A/B, speed up by clustering and designing separately
        if [[ $dataset == influenza_a_* ]] || [[ $dataset == influenza_b_* ]]; then
            cds="--cluster-and-design-separately 0.1"
            is_influenza_ab=true
        else
            cds=""
            is_influenza_ab=false
        fi

        for m in "${MISMATCHES_TO_TRY[@]}"; do
            # --filter-with-lsh-hamming should be commensurate with, but not
            # greater than, m. Use m
            filter_with_lsh_hamming=$m

            for e in "${EXTENSIONS_TO_TRY[@]}"; do
                for i in "${ISLAND_OF_EXACT_MATCHES_TO_TRY[@]}"; do
                    if [ "$is_influenza_ab" = true ]; then
                        # For speed, only compute for large values
                        if [ "$m" -lt "4" ]; then
                            continue
                        fi
                        if [ "$e" -lt "25" ]; then
                            continue
                        fi
                    fi

                    outdir="dataset-runs/$dataset/m${m}-e${e}-i${i}"
                    mkdir -p $outdir
                    outlogfn="$outdir/log.err"
                    outprobesfn="$outdir/probes.fasta"
                    outanalysisfn="$outdir/coverage-analysis.tsv"
                    taxidaccdir="$outdir/taxid-acc"
                    mkdir -p $taxidaccdir

                    if ! [[ -f "${outprobesfn}.gz" && -s "${outprobesfn}.gz" ]]; then
                        # Create a command that writes the number of probes to stdout,
                        # Also, toss the stderr (log) when the job is done
                        cmd="design.py $dataset_input -pl 75 -ps 25 -l 75 -m $m -e $e --island-of-exact-match $i --filter-with-lsh-hamming $filter_with_lsh_hamming --expand-n 0 $cds -o $outprobesfn --write-analysis-to-tsv $outanalysisfn --write-taxid-acc $taxidaccdir --max-num-processes 8 --ncbi-api-key $NCBI_API_KEY --verbose &> $outlogfn; rm $outlogfn; gzip $outprobesfn; gzip $outanalysisfn"
                        echo "$cmd" >> $COMMANDS
                    fi
                done
            done
        done
    done

elif [[ $1 == "param-exploration-parallel-run" ]]; then
    # Run commands that run CATCH for determining number of probes
    # at different choices of parameters
    # Input: $COMMANDS file listing commands from `param-exploration-make-commands`
    # Output: Result of calling those commands - i.e., number of probes written to
    #         dataset-runs/[dataset]/[param choice]/probes.fasta.gz
    # and analysis, etc. for each param choice

    if [ -z "$2" ]; then
        echo "FATAL: Unknown number of commands to run in parallel"
        exit 1
    fi

    NJOBS="$2"

    conda activate /ebs/hybsel-design-runs/tools/envs/catch-prod

    parallel --jobs $NJOBS --no-notice --progress < $COMMANDS

    conda deactivate

elif [[ $1 == "param-exploration-combine-counts" ]]; then
    # Combine probe counts to construct a table of the number of
    # probes for each choice of parameters for each dataset
    # Input: Probe counts given by dataset-runs/[dataset]/[param choice]/probes.fasta.gz
    # Output: Table combining all the probe counts, across datasets and
    #         choices of parameters - written to pooled/probe-counts.tsv
    mkdir -p pooled
    outfn="pooled/probe-counts.tsv"

    # Print a header
    echo -e "dataset\tmismatches\tcover_extension\tisland_of_exact_match\tnum_probes" > $outfn

    # Print the number of probes for each combination of parameters
    for dataset in "${DATASETS[@]}"; do
        for m in "${MISMATCHES_TO_TRY[@]}"; do
            for e in "${EXTENSIONS_TO_TRY[@]}"; do
                for i in "${ISLAND_OF_EXACT_MATCHES_TO_TRY[@]}"; do
                    outprobesfn="dataset-runs/$dataset/m${m}-e${e}-i${i}/probes.fasta.gz"
                    # Only print if the file exists
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
    num_probes=$2
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

    # Run pool.py 100 times (each has a random initial start), and find the
    # overall minimum loss; use the parameters for this run
    conda activate /ebs/hybsel-design-runs/tools/envs/catch-prod
    mkdir -p /tmp/opt-pooled-params
    min_loss="999999999.0"
    min_i="0"
    for i in $(seq 1 100); do
        pool.py pooled/probe-counts.2params.tsv $num_probes /tmp/opt-pooled-params/${i}.params --round-params 1 25 --loss-coeffs 1 0.01 &> /tmp/opt-pooled-params/${i}.out
        loss=$(cat /tmp/opt-pooled-params/${i}.out | grep '^Loss:' | awk '{print $2}')
        if (( $(echo "$loss < $min_loss" | bc -l ) )); then
            min_loss="$loss"
            min_i="$i"
        fi
    done
    cp /tmp/opt-pooled-params/${min_i}.params pooled/design-${num_probes}/params.tsv
    cp /tmp/opt-pooled-params/${min_i}.out pooled/design-${num_probes}/pool.out
    rm -rf /tmp/opt-pooled-params
    conda deactivate

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
        probesfn="dataset-runs/$dataset/m${m}-e${e}-i${i}/probes.fasta.gz"
        zcat $probesfn >> $pooledprobesfn
    done < <(cat $paramsfn | tail -n +2)
    gzip $pooledprobesfn

    # Combine all the coverage analyses
    coverageanalysisfn="pooled/design-${num_probes}/coverage-analysis.tsv"
    while read line; do
        dataset=$(echo "$line" | awk '{print $1}')
        m=$(echo "$line" | awk '{print $2}')
        e=$(echo "$line" | awk '{print $3}')
        i=$(echo "$line" | awk '{print $4}')
        # Get the header of the analysis tsv
        analysisfn="dataset-runs/$dataset/m${m}-e${e}-i${i}/coverage-analysis.tsv.gz"
        zcat $analysisfn | head -n 1 > $coverageanalysisfn
    done < <(cat $paramsfn | tail -n +2 | head -n 1)
    # Concatenate all analysis tsvs, leaving out the header
    while read line; do
        dataset=$(echo "$line" | awk '{print $1}')
        m=$(echo "$line" | awk '{print $2}')
        e=$(echo "$line" | awk '{print $3}')
        i=$(echo "$line" | awk '{print $4}')
        analysisfn="dataset-runs/$dataset/m${m}-e${e}-i${i}/coverage-analysis.tsv.gz"
        zcat $analysisfn | tail -n +2 >> $coverageanalysisfn
    done < <(cat $paramsfn | tail -n +2)
    gzip $coverageanalysisfn
else
    echo "FATAL: Unknown subcommand"
    exit 1
fi
