#!/bin/python3
"""For each fasta file in a results directory, generate a bsub command to
run the 'n_expansion_filter'.

Note that if the probe sequences in the input fasta files have adapters,
this probably won't work correctly because the 'filter_from_fasta' will
filter all probes. In that case, the adapters should be removed from the
probe sequences before they are used as input.
"""

import argparse
import os
import re

FASTA_PATTERN = re.compile('mismatches_([0-9]+)-coverextension_([0-9]+)\.fasta')


def fasta_iter(results_dir):
    for dir in os.listdir(results_dir):
        dataset_results_path = os.path.join(results_dir, dir)
        if not os.path.isdir(dataset_results_path):
            # there are some other files in the results path (e.g.,
            # 'datasets.txt'); only look at the directories
            continue

        # dir gives a virus/dataset name
        dataset = dir
        for fn in os.listdir(dataset_results_path):
            # match fasta files; the parameters are part of the
            # file name
            m = FASTA_PATTERN.match(fn)
            if m:
                # fn is a fasta file
                mismatches = int(m.group(1))
                cover_extension = int(m.group(2))
                fn_path = os.path.abspath(os.path.join(dataset_results_path,
                                                       fn))
                path_prefix = fn_path[:-len('.fasta')]
                yield (mismatches, cover_extension, dataset, path_prefix)

def job_completed_successfully(out_path):
    if not os.path.isfile(out_path):
        return False
    with open(out_path) as f:
        for line in f:
            if 'Successfully completed' in line:
                return True
    return False

def main(args):
    for run in fasta_iter(args.results_dir):
        mismatches, cover_extension, dataset, path_prefix = run
        in_fasta = path_prefix + '.fasta'
        out_fasta = path_prefix + '.n_expanded.fasta'
        bsub_out = path_prefix + '.n_expanded.out'
        if os.path.isfile(out_fasta) and job_completed_successfully(bsub_out):
            # already ran, so skip it
            continue

        mem = 16
        cmd = ["bsub"]
        cmd += ["-o", bsub_out]
        cmd += ["-q", "week"]
        cmd += ["-R", "\"rusage[mem=" + str(mem) + "]\""]
        cmd += ["-P", "hybseldesign"]
        cmd += ["python", "bin/make_probes.py"]
        cmd += ["--probe_length", "75"]
        cmd += ["--probe_stride", "25"]
        cmd += ["--mismatches", str(mismatches)]
        cmd += ["--island_of_exact_match", "30"]
        cmd += ["--cover_extension", str(cover_extension)]
        cmd += ["--dataset", dataset]
        cmd += ["--filter_from_fasta", in_fasta]
        cmd += ["--skip_set_cover"]
        cmd += ["--skip_adapters"]
        cmd += ["--skip_reverse_complements"]
        cmd += ["--expand_n"]
        cmd += ["-o", out_fasta]
        cmd += ["--print_analysis"]
        cmd += ["--verbose"]
        print(' '.join(cmd))


if __name__ == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--results_dir', '-i', required=True)
    args = argparse.parse_args()

    main(args)
