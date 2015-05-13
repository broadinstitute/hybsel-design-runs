#!/bin/python3

import argparse
import itertools
from collections import defaultdict
from os import listdir
from os.path import join

import numpy as np

def read_input(args):
    dataset_coverages = {}
    for filename in listdir(args.input):
        # Only read ".analysis.tsv" files
        if not filename.endswith(".analysis.tsv"):
            continue
        dataset = filename.replace('.analysis.tsv', '')
        coverages = []
        with open(join(args.input, filename)) as f:
            for i, line in enumerate(f):
                if i == 0:
                    # Skip header line
                    continue
                ls = line.rstrip().split('\t')
                genome_name = ls[0]
                if genome_name.endswith("(rc)"):
                    # Skip reverse complements
                    continue
                # Take the average coverage/depth over
                # unambiguous bases
                covg = float(ls[5])
                coverages += [covg]
        dataset_coverages[dataset] = coverages
    return dataset_coverages

def min_and_max_covg(dataset_coverages):
    min_covg = min(itertools.chain(*dataset_coverages.values()))
    max_covg = max(itertools.chain(*dataset_coverages.values()))
    # Drop min_covg down to nearest 0.5
    if int(min_covg) + 0.5 <= min_covg:
        min_covg = int(min_covg) + 0.5
    else:
        min_covg = int(min_covg)
    # Raise max_covg up to the nearest 0.5
    if max_covg != int(max_covg):
        if int(max_covg) + 0.5 > max_covg:
            max_covg = int(max_covg) + 0.5
        else:
            max_covg = int(max_covg) + 1
    return min_covg, max_covg

def main(args):
    dataset_coverages = read_input(args)
    min_covg, max_covg = min_and_max_covg(dataset_coverages)
    # Make bins from min_covg to max_covg, inclusive, separated
    # by 0.5
    bins = np.arange(min_covg, max_covg + 0.5, 0.5)

    # Print header (bin cutoffs)
    header = ["Dataset"]
    for i in xrange(len(bins) - 1):
        header += [str(bins[i]) + "-" + str(bins[i+1])]
    print '\t'.join(header)

    # Print histogram for each dataset
    for dataset in sorted(dataset_coverages.keys()):
        a = dataset_coverages[dataset]
        hist, _ = np.histogram(a, bins=bins)
        # values in hist give the number of genomes falling into
        # each bin; convert these to the fraction of genomes
        hist_frac = [float(n) / len(a) for n in hist]

        row = [dataset] + [str(f) for f in hist_frac]
        print '\t'.join(row)


if __name__ == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--input', '-i', required=True,
        help=("Folder containing coverage analyses by dataset "
              "(files of the form dataset.analysis.tsv)"))
    args = argparse.parse_args()

    main(args)
