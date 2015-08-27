#!/bin/python3
"""Analyze the sequence composition of a set of probes.
"""

import argparse
from collections import defaultdict
import operator

import numpy as np
from scipy.stats import binom
from scipy.stats import hypergeom

import plotly.plotly as py
from plotly.graph_objs import *


def read_probe_seqs(args):
    probe_seqs = {}
    with open(args.probe_seqs) as f:
        curr_header = None
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                curr_header = line[1:]
                continue
            if ("reverse complement of" in curr_header and
                    args.skip_reverse_complement_probes):
                continue
            probe_name = curr_header.split(' | ')[0]
            seq = line
            if args.skip_adapters:
                for adapter in args.adapters:
                    if seq.startswith(adapter):
                        seq = seq[len(adapter):]
                    if seq.endswith(adapter):
                        seq = seq[:-len(adapter)]
            probe_seqs[probe_name] = seq
    return probe_seqs


def read_probe_list(fn):
    names = []
    with open(fn) as f:
        for line in f:
            names += [line.rstrip()]
    return names


def get_seqs(probe_seqs, probes_to_compute):
    return [probe_seqs[p] for p in probes_to_compute]


def find_kmer_counts(seqs, k):
    # Count number of each kmer across seqs
    kmer_counts = defaultdict(int)
    for seq in seqs:
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:(i+k)]
            kmer_counts[kmer] += 1

    num_kmers = sum(kmer_counts[kmer] for kmer in kmer_counts.keys())

    return kmer_counts, num_kmers


def find_kmer_probe_counts(seqs, k):
    # Count number of probes with each kmer (i.e., unlike find_kmer_counts(),
    # a kmer is only counted once if it is contained more than once in
    # a single probe)
    kmer_probe_counts = defaultdict(int)
    for seq in seqs:
        kmers = set()
        for i in range(0, len(seq) - k + 1):
            kmers.add(seq[i:(i+k)])
        for kmer in kmers:
            kmer_probe_counts[kmer] += 1

    return kmer_probe_counts


def print_kmer_freqs(seqs, k=1):
    kmer_counts, num_kmers = find_kmer_counts(seqs, k)
    for kmer in sorted(kmer_counts.keys()):
        frac = float(kmer_counts[kmer]) / num_kmers
        print(kmer + ' : ' + "{:.2%}".format(frac))


def print_kmer_probe_freqs(seqs, k=1):
    kmer_probe_counts = find_kmer_probe_counts(seqs, k)
    for kmer in sorted(kmer_probe_counts.keys()):
        frac = float(kmer_probe_counts[kmer]) / len(seqs)
        print(kmer + ' : ' + "{:.2%}".format(frac))


def find_significant_kmers(foreground_seqs, background_seqs, k=1,
        alpha=0.05, method='hypergeom'):
    # If method is 'hypergeom':
    # Evaluate significance using a hypergeometric test where kmers from the
    # background sequences make up the total bin, and the drawn objects are
    # kmers in the foreground sequences
    #
    # If method is 'binomial':
    # Evaluate significance using a binomial test where each experiment is
    # pulling out a kmer and a success for the experiment is if the
    # kmer is a desired kmer

    foreground_kmer_counts, foreground_num_kmers = \
        find_kmer_counts(foreground_seqs, k)
    background_kmer_counts, background_num_kmers = \
        find_kmer_counts(background_seqs, k)
    all_kmers = set(foreground_kmer_counts.keys()).union(
        set(background_kmer_counts.keys()))

    kmer_pval = {}
    for kmer in all_kmers:
        num_kmer_success_draws = foreground_kmer_counts[kmer]
        num_draws = foreground_num_kmers
        total_kmer_successes = background_kmer_counts[kmer]
        total_num_objects = background_num_kmers
        
        # Survival function (sf) is (1-cdf)
        if method == 'hypergeom':
            pval = hypergeom.sf(num_kmer_success_draws - 1,
                                total_num_objects,
                                total_kmer_successes,
                                num_draws)
        elif method == 'binomial':
            success_prob = float(total_kmer_successes) / total_num_objects
            pval = binom.sf(num_kmer_success_draws - 1,
                            num_draws,
                            success_prob)
        else:
            raise ValueError("Unknown method")

        kmer_pval[kmer] = pval

    # Perform Bonferroni correction
    kmer_pval_corrected = {kmer: min(1.0, kmer_pval[kmer] * len(all_kmers))
                           for kmer in kmer_pval.keys()}

    # Sort by p-value
    sorted_kmer = sorted(kmer_pval_corrected.items(), key=operator.itemgetter(1))
    significant_kmer = [(kmer, pval) for kmer, pval in sorted_kmer if pval < alpha]
    return significant_kmer


def find_significant_kmers_by_probe(foreground_seqs, background_seqs, k=1,
        alpha=0.05):
    # Evaluate significance using a hypergeometric test where the background
    # sequences make up the total bin, and the drawn objects are foreground
    # sequences (which may or may not contain a particular kmer, determining
    # whether the draw is a 'success')

    foreground_kmer_probe_counts = find_kmer_probe_counts(foreground_seqs, k)
    background_kmer_probe_counts = find_kmer_probe_counts(background_seqs, k)
    all_kmers = set(foreground_kmer_probe_counts.keys()).union(
        set(background_kmer_probe_counts.keys()))

    kmer_pval = {}
    for kmer in all_kmers:
        num_kmer_success_draws = foreground_kmer_probe_counts[kmer]
        num_draws = len(foreground_seqs)
        total_kmer_successes = background_kmer_probe_counts[kmer]
        total_num_objects = len(background_seqs)

        if num_kmer_success_draws < 0.001*num_draws:
            # heuristic to skip this kmer (hypergeom.sf is slow)
            continue

        # Survival function (sf) is (1-cdf)
        pval = hypergeom.sf(num_kmer_success_draws - 1,
                            total_num_objects,
                            total_kmer_successes,
                            num_draws)

        kmer_pval[kmer] = pval

    # Perform Bonferroni correction
    kmer_pval_corrected = {kmer: min(1.0, kmer_pval[kmer] * len(all_kmers))
                           for kmer in kmer_pval.keys()}

    # Sort by p-value
    sorted_kmer = sorted(kmer_pval_corrected.items(), key=operator.itemgetter(1))
    significant_kmers = [(kmer, pval) for kmer, pval in sorted_kmer if pval < alpha]

    # Add some more info to each kmer
    for i in range(len(significant_kmers)):
        kmer, pval = significant_kmers[i]
        num_foreground_with_kmer = sum(1 for seq in foreground_seqs if kmer in seq)

        # Compute the fraction of foreground sequences with this kmer
        frac_of_foreground_with_kmer = num_foreground_with_kmer / float(len(foreground_seqs))

        # Compute the fraction of foreground sequences with this kmer that
        # may 'hairpin' due to this kmer
        rc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        kmer_rc = ''.join([rc.get(b, b) for b in kmer][::-1])
        num_foreground_with_kmer_hairpin = sum(1 for seq in foreground_seqs
            if kmer in seq and kmer_rc in seq)
        frac_of_foreground_with_hairpin = num_foreground_with_kmer_hairpin / \
            float(num_foreground_with_kmer)

        significant_kmers[i] = (kmer, pval, frac_of_foreground_with_kmer,
                                frac_of_foreground_with_hairpin)

    return significant_kmers


def print_significant_kmers(significant_kmers, adapters):
    for kmer_info in significant_kmers:
        kmer, pval, other = kmer_info[0], kmer_info[1], kmer_info[2:]
        print('  ', kmer, '(p = ' + str(pval) + ')', end='')
        adapters_with_probe = [a for a in adapters if kmer in a]
        if adapters_with_probe:
            print(' (in adapters ' + ','.join(adapters_with_probe) + ')', end='')
        print(' [' + str(other) + ']')


def create_gc_content_hist(seqs, name, bin_size=0.02):
    gc_fracs = [float(seq.count('G') + seq.count('C'))/len(seq) for
                seq in seqs]
    bins = np.arange(0, 1 + bin_size, bin_size)
    hist, _ = np.histogram(gc_fracs, bins=bins)
    # Values in hist give number of probes (frequency) with gc_frac
    # in each bin; normalize for the total number of probes
    hist_frac = [float(n) / len(gc_fracs) for n in hist]

    scatter = Scatter(
        x=bins,
        y=hist_frac,
        mode='markers+lines',
        name=name
    )
    return scatter
    

def main(args):
    probe_seqs = read_probe_seqs(args)

    gc_content_data = []

    print('BACKGROUND BASE COMPOSITION')
    print('===========================')
    background_seqs = get_seqs(probe_seqs, probe_seqs.keys())
    print_kmer_freqs(background_seqs)
    gc_content_data += [create_gc_content_hist(background_seqs,
                                               'background')]

    if args.probe_names_foreground:
        print('')
        print('FOREGROUND BASE COMPOSITION')
        print('===========================')

        for f in args.probe_names_foreground:
            print('')
            print(f)
            print('-'*len(f))
            probe_names = read_probe_list(f)
            f_seqs = get_seqs(probe_seqs, probe_names)
            print_kmer_freqs(f_seqs)
            print('Enriched kmers:')
            print_significant_kmers(
                find_significant_kmers_by_probe(f_seqs, background_seqs, k=12),
                args.adapters)
            gc_content_data += [create_gc_content_hist(f_seqs,
                                                       'foreground')]

    if args.plot_gc_content_hist:
        gc_content_data = Data(gc_content_data)
        gc_content_layout = Layout(title='GC Content')
        gc_content_fig = Figure(data=gc_content_data,
                                layout=gc_content_layout)
        py.plot(gc_content_fig, filename='gc_content')


if __name__ == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--probe_seqs', required=True,
        help="FASTA file giving all probe sequences")
    argparse.add_argument('--probe_names_foreground', nargs='+',
        help=("List of txt files, each of which lists names of probes (from "
              "probe_seqs file) to consider as a foreground set"))
    argparse.add_argument('--skip_reverse_complement_probes',
        dest='skip_reverse_complement_probes', action='store_true',
        help=("When set, skip probes in probe_seqs whose name contains '"
              "reverse complement'"))
    argparse.add_argument('--plot_gc_content_hist',
        dest='plot_gc_content_hist', action='store_true',
        help=("When set, use plotly to plot the histogram of GC content "
              "in foregrounds vs. background"))
    argparse.add_argument('--adapters', nargs='+',
        default=['ATACGCCATGCTGGGTCTCC', 'CGTACTTGGGAGTCGGCCAT',
                 'AGGCCCTGGCTGCTGATATG', 'GACCTTTTGGGACAGCGGTG'],
        help=("For use in identifying when a kmer is from an adapter"))
    argparse.add_argument('--skip_adapters',
        dest='skip_adapters', action='store_true',
        help=("When set, remove adapters from probes"))
    args = argparse.parse_args()

    main(args)
