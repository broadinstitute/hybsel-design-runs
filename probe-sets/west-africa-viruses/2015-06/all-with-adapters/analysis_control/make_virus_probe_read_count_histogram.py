#!/bin/python3

import argparse
import itertools
from collections import defaultdict
from os import listdir
from os.path import join

import numpy as np

import plotly.plotly as py
from plotly.graph_objs import *


def read_input(args):
    dataset_read_counts = defaultdict(list)
    with open(args.input) as f:
        for line in f:
            ls = line.rstrip().split('\t')
            probe_name = ls[0]
            dataset = ls[1]
            read_count = int(ls[2])
            dataset_read_counts[dataset].append(read_count)
    return dict(dataset_read_counts)


def plot_hist(read_counts, bins, name):
    hist, _ = np.histogram(read_counts, bins=bins)
    # Values in hist give the number of probes (frequency) for each bin of
    # read counts; normalize for the number of probes for the dataset ((len(read_counts))
    hist_frac = [float(n) / len(read_counts) for n in hist]

    up_to = 0
    cum_frac = 0
    for i in xrange(len(hist_frac)):
        cum_frac += hist_frac[i]
        if cum_frac > args.max_cum_frac:
            up_to = i
            break
    if up_to == 0:
        up_to = len(hist_frac) - 1
    hist = hist[:up_to]
    hist_frac = hist_frac[:up_to]

    x_vals = bins[:up_to]
    y_vals = hist

    print name
    data = Data([
        Bar(x=x_vals, y=y_vals)
    ])
    layout = Layout(
        title=name,
        xaxis=XAxis(title='Reads per probe'),
        yaxis=YAxis(title='Frequency')
    )
    fig = Figure(data=data, layout=layout)
    py.plot(fig, filename=name)


def main(args):
    dataset_read_counts = read_input(args)

    bins = np.arange(0, 10000 + args.bin_size, args.bin_size)

    # Plot histogram for each dataset
    for dataset in sorted(dataset_read_counts.keys()):
        read_counts = dataset_read_counts[dataset]
        plot_hist(read_counts, bins, dataset)

    # Plot histogram with all read counts
    all_read_counts = list(itertools.chain(*dataset_read_counts.values()))
    plot_hist(all_read_counts, bins, 'all')


if __name__ == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--input', '-i', required=True,
        help="File giving probe name, dataset, read count for probe on each line")
    argparse.add_argument('--max_cum_frac', required=True, type=float,
        help="Show bins such that the cumulative fraction of read counts reaches this")
    argparse.add_argument('--bin_size', required=True, type=int,
        help="Width of each bin")
    args = argparse.parse_args()

    main(args)
