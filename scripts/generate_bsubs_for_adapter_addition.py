import os

DATASETS = [
            "chikungunya",
#            "crimean_congo",
#            "dengue",
            (["ebola_zaire", "ebola2014"], "ebola_zaire-with-2014"),
#            "ebola_nonzaire",
#            "gbv_c",
#            "hepatitis_a",
#            "hepatitis_c",
#            "hiv1_without_ltr",
#            "hiv2_without_ltr",
#            "influenza",
            "lassa",
#            "marburg",
#            "measles",
#            "mers",
#            "rhabdovirus",
#            "rift_valley_fever",
#            "sars",
#            "yellow_fever",
           ]

ORIG_RESULTS_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015/"
ADAPTER_RESULTS_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015_limited_with-adapters/"
PARAMS_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015_limited_with-adapters/params.12000.txt"

# To use custom adapter sequences, set ADAPTER_SEQUENCES_PATH to be a file where
# each line is of the following format:
# [dataset name]	[A adapter 5' end]..[A adapter 3' end]	[B adapter 5' end]..[B adapter 3' end]
# (to use default adapter sequences, set ADAPTER_SEQUENCES_PATH to None)
ADAPTER_SEQUENCES_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015_limited_with-adapters/adapters.txt"

PARAMS = {}
with open(PARAMS_PATH) as f:
    for line in f:
        ls = line.rstrip().split('\t')
        dataset = ls[0]
        params = eval(ls[1])
        PARAMS[dataset] = params

if ADAPTER_SEQUENCES_PATH != None:
    # Use custom adapters
    ADAPTER_SEQUENCES = {}
    with open(ADAPTER_SEQUENCES_PATH) as f:
        for line in f:
            ls = line.rstrip().split('\t')
            dataset = ls[0]
            adapter_a = tuple(ls[1].split('..'))
            adapter_b = tuple(ls[2].split('..'))
            ADAPTER_SEQUENCES[dataset] = (adapter_a, adapter_b)
else:
    ADAPTER_SEQUENCES = None

for dataset in DATASETS:
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

    mismatches, cover_extension = PARAMS[name]
    orig_fasta_path = (ORIG_RESULTS_PATH + name + "/mismatches_" + str(mismatches) +
                       "-coverextension_" + str(cover_extension) + ".fasta")
    path = ADAPTER_RESULTS_PATH + name

    mem = 4
    cmd = ["bsub"]
    cmd += ["-o", path +  ".out"]
    cmd += ["-q", "week"]
    cmd += ["-R", "\"rusage[mem=" + str(mem) + "]\""]
    cmd += ["-P", "hybseldesign"]
    cmd += ["python", "bin/make_probes.py"]
    cmd += ["--mismatches", str(mismatches)]
    cmd += ["--lcf_thres", "100"]
    cmd += ["--cover_extension", str(cover_extension)]
    cmd += ["--dataset"] + dataset
    cmd += ["--filter_from_fasta", orig_fasta_path]
    cmd += ["--skip_set_cover"]
    if ADAPTER_SEQUENCES != None:
        # Use custom adapters
        adapter_a, adapter_b = ADAPTER_SEQUENCES[name]
        cmd += ["--adapter_a"] + list(adapter_a)
        cmd += ["--adapter_b"] + list(adapter_b)
    cmd += ["--print_analysis"]
    cmd += ["-o", path + ".fasta"]
    cmd += ["--verbose"]
    print ' '.join(cmd)
