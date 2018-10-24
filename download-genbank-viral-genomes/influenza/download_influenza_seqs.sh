#!/bin/bash

# Run from ../
python -u download_dataset_fastas.py -dl influenza/influenza.datasets.txt -gl influenza/influenza.acc-list.txt --consolidate-segmented-genomes-into-one-fasta --gzip-fastas -o all

