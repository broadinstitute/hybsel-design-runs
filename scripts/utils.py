#!/bin/python3

import os
import re
import subprocess

FASTA_PATTERN = re.compile('mismatches_([0-9]+)-coverextension_([0-9]+)\.fasta')


def count_probes(fn):
    # return the number of probes in the fasta file fn
    cmd = ["grep", "'>'", fn, "|", "wc", "-l"]
    out = subprocess.check_output(' '.join(cmd), shell=True)
    if not isinstance(out, str):
        out = out.decode('utf-8')
    out = out.rstrip()
    return int(out)


def read_probe_counts(args, skip=["datasets.txt"]):
    probe_counts = {}
    for dir in os.listdir(args.results_dir):
        if os.path.isfile(dir) and dir in skip:
            continue
        else:
            assert os.path.isdir(dir)
        # dir gives a virus/dataset name
        dataset = dir
        dataset_results_path = os.path.join(args.results_dir, dir)
        d = {}
        for fn in os.listdir(dataset_results_path):
            # match fasta files; the parameters are part of the
            # file name
            m = FASTA_PATTERN.match(fn)
            if m:
                # fn is a fasta file
                mismatches = int(m.group(1))
                cover_extension = int(m.group(2))
                fn_path = os.path.join(dataset_results_path, fn)
                probe_count = count_probes(fn_path)
                d[(mismatches, cover_extension)] = probe_count
        probe_counts[dataset] = d
    return probe_counts

