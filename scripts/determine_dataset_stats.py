#!/bin/python3
"""Output stats associated with each hybsel_design dataset.
"""

import argparse
import os

import seq_io


def compute_unsegmented_dataset_stats(data_dir, fasta_fn):
    fasta = seq_io.read_fasta(os.path.join(data_dir, fasta_fn))

    num_genomes = len(fasta)
    num_seqs = num_genomes
    avg_seq_len = sum([len(seq) for seq in fasta.values()]) / float(num_seqs)

    return (num_genomes, num_seqs, avg_seq_len)


def compute_segmented_dataset_stats(data_dir, dataset_folder):
    num_genomes = 0
    num_seqs = 0
    total_seq_len = 0

    for fasta_fn in os.listdir(os.path.join(data_dir, dataset_folder)):
        fasta_fn_path = os.path.join(data_dir, dataset_folder, fasta_fn)
        assert os.path.isfile(fasta_fn_path)

        fasta = seq_io.read_fasta(fasta_fn_path)
        
        num_genomes += 1
        num_seqs += len(fasta)
        total_seq_len += sum([len(seq) for seq in fasta.values()])

    avg_seq_len = total_seq_len / float(num_seqs)

    return (num_genomes, num_seqs, avg_seq_len)


def compute_dataset_stats(data_dir):
    stats = {}
    for fn in os.listdir(data_dir):
        fn_path = os.path.join(data_dir, fn)
        if os.path.isfile(fn_path):
            assert fn.endswith('.fasta')
            dataset_name = fn[:-len('.fasta')]
            stats[dataset_name] = compute_unsegmented_dataset_stats(data_dir,
                                                                    fn)
        else:
            assert os.path.isdir(fn_path)
            dataset_name = fn
            stats[dataset_name] = compute_segmented_dataset_stats(data_dir,
                                                                  fn)
    return stats


def main(args):
    stats = compute_dataset_stats(args.datasets_data_dir)
    for dataset in sorted(stats.keys()):
        print('\t'.join([dataset] + [str(x) for x in stats[dataset]]))


if __name__ == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--datasets_data_dir', required=True,
        help="Folder with hybsel_design datasets data")
    args = argparse.parse_args()

    main(args)
