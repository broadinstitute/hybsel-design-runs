import os

DATASETS = [
            "chikungunya",
            "crimean_congo",
            "dengue",
            (["ebola_zaire", "ebola2014"], "ebola_zaire-with-2014"),
            "ebola_nonzaire",
            "gbv_c",
            "hepatitis_a",
            "hepatitis_c",
            "hiv1_without_ltr",
            "hiv2_without_ltr",
            "influenza",
            "lassa",
            "marburg",
            "measles",
            "mers",
            "rhabdovirus",
            "rift_valley_fever",
            "sars",
            "yellow_fever",
           ]

ORIG_RESULTS_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015/"
ADAPTER_RESULTS_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015_with-adapters/"
PARAMS_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015_with-adapters/params.90000.txt"

PARAMS = {}
with open(PARAMS_PATH) as f:
    for line in f:
        ls = line.rstrip().split('\t')
        dataset = ls[0]
        params = eval(ls[1])
        PARAMS[dataset] = params

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
    cmd += ["--print_analysis"]
    cmd += ["-o", path + ".fasta"]
    cmd += ["--verbose"]
    print ' '.join(cmd)
