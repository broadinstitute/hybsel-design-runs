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

RESULTS_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015/"

# parameter space is (mismatches, cover_extension)
PARAMETER_SPACE = [(mismatches, cover_extension)
                   for mismatches in range(0, 10)
                   for cover_extension in range(0, 51, 10)]

def mem_requested(dataset_name, mismatches, cover_extension):
    big_datasets = ["hepatitis_c", "influenza", "hiv1_without_ltr"]
    if dataset_name in big_datasets:
        if mismatches >= 7:
            return 24
        elif mismatches >= 4:
            return 16
        else:
            return 8
    else:
        if mismatches >= 7:
            return 16
        elif mismatches >= 4:
            return 8
        else:
            return 4

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
    # Make the directory for this dataset's results
    if not os.path.exists(RESULTS_PATH + name):
        os.makedirs(RESULTS_PATH + name)
    for params in PARAMETER_SPACE:
        mismatches, cover_extension = params
        path = (RESULTS_PATH + name + "/mismatches_" + str(mismatches) +
                "-coverextension_" + str(cover_extension))
        mem = mem_requested(name, mismatches, cover_extension)
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
        cmd += ["--skip_adapters"]
        cmd += ["--print_analysis"]
        cmd += ["--write_analysis_to_tsv", path + ".analysis.tsv"]
        cmd += ["--write_sliding_window_coverage", path + ".covg"]
        cmd += ["-o", path + ".fasta"]
        cmd += ["--verbose"]
        print ' '.join(cmd)
