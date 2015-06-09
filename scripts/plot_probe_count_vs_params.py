#!/bin/python3

import argparse
from collections import defaultdict

import plotly.plotly as py
from plotly.graph_objs import *

import utils


def divide_by_cover_extension(kv):
    kv_by_cover_extension = defaultdict(list)
    for param, probe_count in kv.items():
        mismatches, cover_extension = param
        kv_by_cover_extension[cover_extension].append((mismatches, probe_count))
    return dict(kv_by_cover_extension)


def plot_dataset(dataset, kv_by_cover_extension):
    traces = []
    for cover_extension in sorted(kv_by_cover_extension.keys()):
        vals_sorted = sorted(kv_by_cover_extension[cover_extension])
        # mismatches are x values
        x = [v[0] for v in vals_sorted]
        # probe counts are y values
        y = [v[1] for v in vals_sorted]
        trace = Scatter(
            x=x,
            y=y,
            mode='lines+markers',
            name=str(cover_extension),
            line=Line(
                shape='linear'
            )
        )
        traces += [trace]

    data = Data(traces)
    layout = Layout(
        legend=Legend(
            y=0.5,
            font=Font(
                size=16
            ),
            yref='paper'
        ),
        title=dataset,
        xaxis=XAxis(
            title="Mismatches"
        ),
        yaxis=YAxis(
            title="Probe count"
        )
    )
    fig = Figure(data=data, layout=layout)
    py.image.save_as(fig, os.path.join(args.output_dir, dataset+".pdf"))


def determine_mismatch_count(args, kv_by_cover_extension):
    mismatch_count = {}
    for cover_extension in sorted(kv_by_cover_extension.keys()):
        vals_sorted = sorted(kv_by_cover_extension[cover_extension])
        mismatch_pick = None
        if args.mismatch_pick_method == "percent_change_against_next":
            # pick mismatch parameter by looking at the percent change in
            # probe count for mismatch M versus M+1, and pick M if that
            # change is too small
            for i in xrange(len(vals_sorted)):
                mismatches, probe_count = vals_sorted[i]
                if i == len(vals_sorted) - 1:
                    # no more to consider
                    mismatch_pick = mismatches
                    break
                else:
                    next_probe_count = vals_sorted[i+1][1]
                    frac_change = (next_probe_count - probe_count) / float(probe_count)
                    if frac_change > args.mismatch_change_cutoff:
                        # note that frac_change is negative (denoting a decrease in
                        # probe count -- i.e., next_probe_count should be < than
                        # probe_count) and args.mismatch_change_cutoff is negative.
                        # the change to the next probe count is  *too little* (does
                        # not decrease by enough), so use this as the number of
                        # mismatches
                        mismatch_pick = mismatches
                        break
        elif args.mismatch_pick_method == "percent_of_original_against_next":
            # pick mismatch parameter M if the following holds: let P be the
            # number of probes that M yields divided by the number of probes with
            # 0 mismatches (original); compute P with M+1 mismatches and call that
            # P'; if (P - P') is too little, use M as the number of mismatches
            # this relates more to what we visualize looking at the plots
            orig_probe_count = vals_sorted[0][1]
            for i in xrange(len(vals_sorted)):
                mismatches, probe_count = vals_sorted[i]
                if i == len(vals_sorted) - 1:
                    # no more to consider
                    mismatch_pick = mismatches
                    break
                else:
                    probe_count_prct = float(probe_count) / orig_probe_count
                    next_probe_count = vals_sorted[i+1][1]
                    next_probe_count_prct = float(next_probe_count) / orig_probe_count
                    prct_change = next_probe_count_prct - probe_count_prct
                    if prct_change > args.mismatch_change_cutoff:
                        # note that prct_change is negative (denoting a decrease in
                        # probe count -- i.e., next_probe_count should be < than
                        # probe_count) and args.mismatch_change_cutoff is negative.
                        # the change to the next probe count is  *too little* (does
                        # not decrease by enough), so use this as the number of
                        # mismatches
                        mismatch_pick = mismatches
                        break
        else:
            raise ValueError("Unknown mismatch_pick_method")
        mismatch_count[cover_extension] = mismatch_pick
    return mismatch_count


def main(args):
    probe_counts = utils.read_probe_counts(args)

    for dataset in sorted(probe_counts.keys()):
        d = probe_counts[dataset]
        print "On", dataset
        probe_counts_by_cover_extension = divide_by_cover_extension(d)
        if args.plot:
            print "  Plotting"
            plot_dataset(dataset, probe_counts_by_cover_extension)
        if args.show_mismatch_picks:
            print ' ', determine_mismatch_count(args, probe_counts_by_cover_extension)


if __name__ == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--results_dir', '-i', required=True)
    argparse.add_argument('--output_dir', '-o',
        default='/broad/hptmp/hmetsky/')
    argparse.add_argument('--mismatch_change_cutoff', '-t', type=float,
        default=-0.05)
    argparse.add_argument('--mismatch_pick_method', '-m',
        default="percent_of_original_against_next")
    argparse.add_argument('--plot', dest='plot', action='store_true')
    argparse.add_argument('--show_mismatch_picks',
        dest='show_mismatch_picks', action='store_true')
    args = argparse.parse_args()

    main(args)
