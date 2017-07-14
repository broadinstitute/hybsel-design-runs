import os

DATASETS = [
            "chikungunya",
            "crimean_congo",
            "dengue",
            "ebola_zaire-with-2014",
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

FASTA_RESULT_PATH = "/home/unix/hmetsky/viral/viral-work/results/hybsel_design/viral-probe-set_06-2015_with-adapters/"

SKIP_REVERSE_COMPLEMENT_PROBES = True

for dataset in DATASETS:
    fasta_path = os.path.join(FASTA_RESULT_PATH, dataset + ".fasta")
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                # Header
                if "reverse complement of" in line and SKIP_REVERSE_COMPLEMENT_PROBES:
                    continue
                probe_name = line[1:].split(' | ')[0]
                print(probe_name + "\t" + dataset)
