import argparse
import os
from subprocess import call
import sys

import numpy as np

ANALYZE_BIN = "bin/analyze_probe_coverage.py"

DEFAULT_TMP_DIR = "/home/unix/hmetsky/tmp/analyzeprobes/"

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

# parameter space is (mismatches, lcf_thres, cover_extension)
PARAMETER_SPACE = [(mismatches, lcf_thres, cover_extension)
                   for mismatches in range(0, 17)
                   for lcf_thres in [100]
                   for cover_extension in range(0, 101, 10)]


def analysis_path(tmp_dir, dataset_name, mismatches,
                  lcf_thres, cover_extension):
    path = os.path.join(tmp_dir,
                        (dataset_name + "." + 
                         "mismatches_" + str(mismatches) + "." + 
                         "lcfthres_" + str(lcf_thres) + "." +
                         "coverextension_" + str(cover_extension)))
    return path


def iter_dataset():
    # yield (dataset, name)
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
        yield (dataset, name)


def iter_dataset_and_params():
    # yield (dataset, name, params)
    for dataset, name in iter_dataset():
        for params in PARAMETER_SPACE:
            yield (dataset, name, params)


def run_analysis():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp_dir", default=DEFAULT_TMP_DIR,
                        help="tmp directory")
    parser.add_argument("-f", "--probes_fasta", required=True,
                        help="path to fasta with probes")
    args = parser.parse_args(sys.argv[2:])

    for dataset, name, params in iter_dataset_and_params():
        mismatches, lcf_thres, cover_extension = params
        path = analysis_path(args.tmp_dir, name, mismatches, lcf_thres,
                             cover_extension)
        path_analysis = path + ".analysis.tsv"
        if os.path.exists(path_analysis):
            # skip b/c the output file already exists
            continue

        num_processes = 4
        cmd = ["bsub"]
        cmd += ["-o", path +  ".out"]
        cmd += ["-q", "forest"]
        cmd += ["-P", "hybseldesign"]
        cmd += ["-n", str(num_processes)]
        cmd += ["-R", "\"span[hosts=1]\""]
        cmd += ["python", ANALYZE_BIN]
        cmd += ["--mismatches", str(mismatches)]
        cmd += ["--lcf_thres", "100"]
        cmd += ["--cover_extension", str(cover_extension)]
        cmd += ["--dataset"] + dataset
        cmd += ["--probes_fasta", args.probes_fasta]
        cmd += ["--write_analysis_to_tsv", path_analysis]
        cmd += ["--max_num_processes", str(num_processes)]
        #call(cmd)
        print(' '.join(cmd))


def clean():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp_dir", default=DEFAULT_TMP_DIR)
    args = parser.parse_args(sys.argv[2:])

    # remove all files (but not subdirs) in the tmp dir
    for f in os.listdir(args.tmp_dir):
        fp = os.path.join(args.tmp_dir, f)
        if os.path.isfile(fp):
            os.unlink(fp)


def read_frac_of_covered_genome(fn):
    # Returns list of fracs of unambig bases covered, one per genome
    # (skips reverse complement genomes)
    fracs = []
    with open(fn) as f:
        next(f) # skip header
        for line in f:
            ls = line.rstrip().split('\t')
            genome_name = ls[0]
            if genome_name.endswith('(rc)'):
                # skip reverse complement genome
                continue
            frac_unambig_covered = float(ls[3])
            fracs += [frac_unambig_covered]
    return fracs


def coverage_is_acceptable(fracs_unambig_covered):
    # input: for each genome, fraction of unambiguous bases that are covered
    # acceptable if the 5'th percentile (i.e., 5% of least covered genomes)
    #  has this frac as > 0.99 (i.e., > 99% covered)
    #  (in other words, 95% of genomes are > 99% covered)
    prctile = np.percentile(fracs_unambig_covered, 5.0)
    return prctile > 0.99


def params_loss(params):
    mismatches, lcf_thres, cover_extension = params
    return (mismatches**2.0 +
            (100.0 - lcf_thres)**2.0 +
            (cover_extension / 5.0)**2.0)


def best_params_giving_acceptable_coverage(covg_for_params):
    # input: dict mapping params (tuple) to fracs_unambig_covered
    best_acceptable = None
    for params in covg_for_params.keys():
        fracs_unambig_covered = covg_for_params[params]
        if coverage_is_acceptable(fracs_unambig_covered):
            if best_acceptable is None:
                best_acceptable = params
            elif params_loss(params) < params_loss(best_acceptable):
                best_acceptable = params
    return best_acceptable


def summarize():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp_dir", default=DEFAULT_TMP_DIR)
    args = parser.parse_args(sys.argv[2:])

    dataset_params_covg = {}
    for dataset, name, params in iter_dataset_and_params():
        mismatches, lcf_thres, cover_extension = params
        path = analysis_path(args.tmp_dir, name, mismatches, lcf_thres,
                             cover_extension)
        path_full = path + ".analysis.tsv"
        if not os.path.exists(path_full):
            print("MISSING FILE", path_full, file=sys.stderr)
            continue
        fracs_unambig_covered = read_frac_of_covered_genome(path_full)

        if name not in dataset_params_covg:
            dataset_params_covg[name] = {}
        dataset_params_covg[name][params] = fracs_unambig_covered

    for dataset_name in sorted(dataset_params_covg.keys()):
        best_acceptable_params = best_params_giving_acceptable_coverage(
            dataset_params_covg[dataset_name])
        print(dataset_name, best_acceptable_params)


def specified_param(percentiles=[0, 5, 25, 50, 75, 95, 100], default_lcf=100):
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp_dir", default=DEFAULT_TMP_DIR)
    parser.add_argument("--param_choices", required=True,
        help=("file giving, on each line, "
              "'[dataset name]  [param choices as tuple]'"))
    args = parser.parse_args(sys.argv[2:])

    dataset_params = {}
    with open(args.param_choices) as f:
        for line in f:
            ls = line.rstrip().split('\t')
            dataset = ls[0]
            params = eval(ls[1])
            dataset_params[dataset] = params

    for dataset_name in sorted(dataset_params.keys()):
        if len(dataset_params[dataset_name]) == 2:
            mismatches, cover_extension = dataset_params[dataset_name]
            lcf_thres = default_lcf
        else:
            mismatches, lcf_thres, cover_extension = dataset_params[dataset_name]
        path = analysis_path(args.tmp_dir, dataset_name, mismatches, lcf_thres,
                             cover_extension)
        path_full = path + ".analysis.tsv"
        if not os.path.exists(path_full):
            print("MISSING FILE", path_full, file=sys.stderr)
            continue
        fracs_unambig_covered = read_frac_of_covered_genome(path_full)

        covgs = np.percentile(fracs_unambig_covered, percentiles)
        covgs = [min(1.0, c) for c in covgs]
        covgs_formatted = [float("{0:.2f}".format(100.0*c)) for c in covgs]
        print(dataset_name, covgs_formatted)


if __name__ == "__main__":
    this_module = sys.modules[__name__]
    
    parser = argparse.ArgumentParser(
        usage="""python analyze_probe_coverage_from_fasta.py <command> [<args>]
        The available commands are:
            clean           delete contents of the tmp dir
            run_analysis    run hybseldesign's coverage analyzer for each
                            dataset and param combination
            summarize       collect information from all of the runs, across
                            all datasets and param combinations
            specified_param summarize coverage information for each dataset,
                            using only a specified parameter choice for
                            each dataset
        """)

    parser.add_argument('command', help="Subcommand to run")
    args = parser.parse_args(sys.argv[1:2])
    if not hasattr(this_module, args.command):
        print("Unknown command")
        parser.print_help()
        exit(1)
    # invoke the function with the name given by 'command'
    getattr(this_module, args.command)()
