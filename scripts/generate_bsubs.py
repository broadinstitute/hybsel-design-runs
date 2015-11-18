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

# if re-running some jobs because they failed due to memory limits
# (submitted with too little memory requested), re-run them with
# double the initial requested amount
DOUBLE_MEM = True

# output of 'bjobs -w | grep RUN'; don't submit commands
# that are running
RUNNING_CMD_LIST = "/home/unix/hmetsky/tmp/running"
RUNNING = []
if RUNNING_CMD_LIST is not None:
    with open(RUNNING_CMD_LIST) as f:
        for line in f:
            RUNNING += [line.rstrip()]


# parameter space is (mismatches, cover_extension)
PARAMETER_SPACE = [(mismatches, cover_extension)
                   for mismatches in range(0, 10)
                   for cover_extension in range(0, 51, 10)]

def mem_requested(num_seqs, avg_seq_len, mismatches):
    cost = num_seqs * avg_seq_len
    if cost > 2 * 10**7:
        if mismatches >= 7:
            mem = 32
        elif mismatches >= 4:
            mem = 24
        else:
            mem = 16
    elif cost > 10**7:
        if mismatches >= 7:
            mem = 24
        elif mismatches >= 4:
            mem = 16
        else:
            mem = 8
    else:
        if mismatches >= 7:
            mem = 8
        elif mismatches >= 4:
            mem = 4
        else:
            mem = 2

    if DOUBLE_MEM:
        mem *= 2
    return mem

def queue_requested(num_seqs, avg_seq_len, mismatches, mem):
    if mem > 64:
        return "forest"
    else:
        cost = num_seqs * avg_seq_len
        if cost > 10**6:
            return "week"
        else:
            return "hour"

def job_completed_successfully(out_path):
    if not os.path.isfile(out_path):
        return False
    with open(out_path) as f:
        for line in f:
            if 'Successfully completed' in line:
                return True
    return False

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

        if job_completed_successfully(path + '.out'):
            continue

        mem = mem_requested(num_seqs, avg_seq_len, mismatches)
        queue = queue_requested(num_seqs, avg_seq_len, mismatches, mem)

        bsub_cmd = ["bsub"]
        bsub_cmd += ["-o", path +  ".out"]
        bsub_cmd += ["-q", queue]
        bsub_cmd += ["-R", "\"rusage[mem=" + str(mem) + "]\""]
        bsub_cmd += ["-P", "hybseldesign"]

        cmd = ["python", "bin/make_probes.py"]
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

        py_cmd = ' '.join(cmd)
        skip_cmd = False
        for r in RUNNING:
            if py_cmd in r:
                skip_cmd = True
                break
        if skip_cmd:
            continue

        print(' '.join(bsub_cmd + cmd))
