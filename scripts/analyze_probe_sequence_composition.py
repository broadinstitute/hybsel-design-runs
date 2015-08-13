#!/bin/python3
"""Analyze the sequence composition of a set of probes.
"""

import argparse

import numpy as np

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
            probe_seqs[probe_name] = line
    return probe_seqs


def read_probe_list(fn):
    names = []
    with open(fn) as f:
        for line in f:
            names += [line.rstrip()]
    return names


def get_seqs(probe_seqs, probes_to_compute):
    return [probe_seqs[p] for p in probes_to_compute]


def print_kmer_composition(seqs, k=1):
    bases = ['A', 'T', 'C', 'G']
    kmers = list(bases)
    for i in xrange(k - 1):
        new_kmers = []
        for kmer in kmers:
            for b in bases:
                new_kmers += [kmer + b]
        kmers = new_kmers

    kmer_counts = {kmer: 0 for kmer in kmers}
    for seq in seqs:
        for kmer in kmer_counts.keys():
            kmer_counts[kmer] += seq.count(kmer)
    num_kmers = sum(kmer_counts[kmer] for kmer in kmer_counts.keys())
    for kmer in sorted(kmer_counts.keys()):
        frac = float(kmer_counts[kmer]) / num_kmers
        print kmer + ' : ' + "{:.2%}".format(frac)


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

    print 'BACKGROUND BASE COMPOSITION'
    print '==========================='
    background_seqs = get_seqs(probe_seqs, probe_seqs.keys())
    print_kmer_composition(background_seqs)
    gc_content_data += [create_gc_content_hist(background_seqs,
                                               'background')]

    if args.probe_names_foreground:
        print ''
        print 'FOREGROUND BASE COMPOSITION'
        print '==========================='

        for f in args.probe_names_foreground:
            print ''
            print f
            print '-'*len(f)
            probe_names = read_probe_list(f)
            f_seqs = get_seqs(probe_seqs, probe_names)
            print_kmer_composition(f_seqs)
            gc_content_data += [create_gc_content_hist(f_seqs,
                                                       'foreground')]

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
    args = argparse.parse_args()

    main(args)
