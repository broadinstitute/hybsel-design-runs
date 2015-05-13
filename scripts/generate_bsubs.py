DATASETS = [
            "chikungunya",
            "dengue",
            (["ebola_zaire", "ebola2014"], "ebola_zaire-with-2014"),
            "ebola_nonzaire",
            "hepatitis_a",
            "hepatitis_c",
            "influenza",
            "lassa",
            "marburg",
            "measles",
            "mers",
            "rhabdovirus",
            "sars",
            "yellow_fever"
           ]

RESULTS_PATH = "~/viral/viral-work/results/hybsel_design/viral-probe-set_05-2015/"

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
    cmd = ["bsub"]
    cmd += ["-o", RESULTS_PATH + name + ".out"]
    cmd += ["-q", "forest"]
    cmd += ["python", "bin/make_probes.py"]
    cmd += ["-m", "3", "-l", "100"]
    cmd += ["--dataset"] + dataset
    cmd += ["--skip_adapters"]
    cmd += ["--print_analysis"]
    cmd += ["--write_analysis_to_tsv", RESULTS_PATH + name + ".analysis.tsv"]
    cmd += ["--write_sliding_window_coverage", RESULTS_PATH + name + ".covg"]
    cmd += ["-o", RESULTS_PATH + name + ".fasta"]
    cmd += ["--verbose"]

    print ' '.join(cmd)
