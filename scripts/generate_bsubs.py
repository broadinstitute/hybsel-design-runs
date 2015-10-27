from collections import OrderedDict
import os


RESULTS_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_all-human-host-viruses/2015-10/"

DATASETS = OrderedDict()
with open(RESULTS_PATH + "datasets.txt") as f:
    for line in f:
        ls = line.split('\t')
        dataset = ls[0]
        stats = (int(ls[1]), int(ls[2]), float(ls[3]))
        DATASETS[dataset] = stats


# parameter space is (mismatches, cover_extension)
PARAMETER_SPACE = [(mismatches, cover_extension)
                   for mismatches in range(0, 10)
                   for cover_extension in range(0, 51, 10)]

def mem_requested(num_seqs, avg_seq_len, mismatches):
    cost = num_seqs * avg_seq_len
    if cost > 10**6:
        if mismatches >= 7:
            return 24
        elif mismatches >= 4:
            return 16
        else:
            return 8
    else:
        if mismatches >= 7:
            return 8
        elif mismatches >= 4:
            return 4
        else:
            return 2

def queue_requested(num_seqs, avg_seq_len, mismatches):
    cost = num_seqs * avg_seq_len
    if cost > 10**6:
        return "forest"
    else:
        return "hour"

for dataset, stats in DATASETS.items():
    if isinstance(dataset, tuple):
        # dataset is a tuple (x, y) where x is an array of datasets
        # and y is the name that should be given to this collection of
        # datasets
        dataset, name = dataset
    else:
        # This is a single dataset: the name should be the same as
        # dataset
        name = dataset
        dataset = [dataset]

    num_genomes, num_seqs, avg_seq_len = stats

    # Make the directory for this dataset's results
    if not os.path.exists(RESULTS_PATH + name):
        os.makedirs(RESULTS_PATH + name)
    for params in PARAMETER_SPACE:
        mismatches, cover_extension = params
        path = (RESULTS_PATH + name + "/mismatches_" + str(mismatches) +
                "-coverextension_" + str(cover_extension))
        mem = mem_requested(num_seqs, avg_seq_len, mismatches)
        queue = queue_requested(num_seqs, avg_seq_len, mismatches)

        cmd = ["bsub"]
        cmd += ["-o", path +  ".out"]
        cmd += ["-q", queue]
        cmd += ["-R", "\"rusage[mem=" + str(mem) + "]\""]
        cmd += ["-P", "hybseldesign"]
        cmd += ["python", "bin/make_probes.py"]
        cmd += ["--mismatches", str(mismatches)]
        cmd += ["--lcf_thres", "100"]
        cmd += ["--cover_extension", str(cover_extension)]
        cmd += ["--dataset"] + dataset
        cmd += ["--skip_adapters"]
        cmd += ["--print_analysis"]
        cmd += ["--write_analysis_to_tsv", path + ".analysis.tsv"]
        cmd += ["--write_sliding_window_coverage", path + ".covg"]
        cmd += ["-o", path + ".fasta"]
        cmd += ["--verbose"]
        print(' '.join(cmd))
