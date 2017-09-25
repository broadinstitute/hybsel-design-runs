"""Remove the LTR regions (from both 5' and 3' ends) from
HIV FASTA sequences.
"""

import argparse
from collections import OrderedDict
from os.path import dirname
from os.path import join

import numpy as np

import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


# 5' and 3' LTR coordinates in HIV-1 from:
#  http://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html
HIV1_LTR5_COORDS = (0, 634)
HIV1_LTR3_COORDS = (9085, 9719)

# 5' and 3' LTR coordinates in HIV-2 from:
#  http://www.hiv.lanl.gov/content/immunology/pdf/2000/intro/GenomeMaps.pdf
HIV2_LTR5_COORDS = (0, 855)
HIV2_LTR3_COORDS = (9504, 10359)


def main(args):
    seqs = seq_io.read_fasta(args.i)

    # Find the reference sequence. For HIV-1, the sequence header
    # contains 'reference'; for HIV-2, the sequence header contains
    # 'isolate BEN'
    if args.s == 'hiv1':
        ref_search_pattern = 'reference'
    elif args.s == 'hiv2':
        ref_search_pattern = 'isolate BEN'
    else:
        raise ValueError("Unknown sequence %s" % args.s)
    num_with_ref = sum(True if ref_search_pattern in name else False
                       for name in seqs.keys())
    assert num_with_ref == 1
    ref_name = [name for name in seqs.keys() if ref_search_pattern in name][0]
    ref_seq = seqs[ref_name]

    # Find the index end of the 5' LTR (i.e., counting gaps)
    if args.s == 'hiv1':
        ltr_end_bp = HIV1_LTR5_COORDS[1]
    elif args.s == 'hiv2':
        ltr_end_bp = HIV2_LTR5_COORDS[1]
    bp = 0
    for i in range(len(ref_seq)):
        if ref_seq[i] != '-':
            if bp == ltr_end_bp - 1:
                ltr_5_end = i + 1
                break
            bp += 1

    # Find the index start of the 3' LTR (i.e., counting gaps)
    if args.s == 'hiv1':
        ltr_start_bp = HIV1_LTR3_COORDS[0]
    elif args.s == 'hiv2':
        ltr_start_bp = HIV2_LTR3_COORDS[0]
    else:
        raise ValueError("Unknown sequence %s" % args.s)
    bp = 0
    for i in range(len(ref_seq)):
        if ref_seq[i] != '-':
            if bp == ltr_start_bp:
                ltr_3_start = i
                break
            bp += 1

    # Write the sequences without LTRs
    seqs_no_ltrs = []
    for name, seq in seqs.items():
        # Ensure the seq (with gaps) has the same length as the
        # reference seq
        assert len(seq) == len(ref_seq)
        # Remove the LTR by chopping seq (which includes the gap
        # character '-') to only include ltr_5_end:ltr_3_start
        seq_no_ltr = seq[ltr_5_end:ltr_3_start]
        # Remove gaps from the sequence
        seq_no_ltr = seq_no_ltr.replace('-', '')
        seqs_no_ltrs += [(name, seq_no_ltr)]
    seq_io.write_fasta(OrderedDict(seqs_no_ltrs), args.o)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True,
        help="Input aligned FASTA file")
    parser.add_argument('-o', required=True,
        help="Output FASTA file")
    parser.add_argument('-s', required=True,
        help="Seq ('hiv1' or 'hiv2')")

    args = parser.parse_args()
    main(args)
