#!/bin/python3

import argparse
import itertools
from collections import defaultdict


def read_input(args):
  kv = defaultdict(dict)
  with open(args.input) as f:
    for line in f:
      ls = line.rstrip().split()
      if len(ls) != 3:
        continue
      x, y, val = int(ls[0]), int(ls[1]), float(ls[2])
      if int(val) == val:
        val = int(val)
      kv[x][y] = val
  return kv


def print_matrix(kv, args):
  row_vals = sorted(list(set(x for x in kv.keys())))
  col_vals = sorted(list(set(itertools.chain(*\
                      [kv[x].keys() for x in kv.keys()]))))

  # print col headers
  print("\t".join(["/"] + [str(x) for x in col_vals]))

  # print rows
  for row in row_vals:
    print("\t".join([str(row)] + \
            [str(x) for x in [kv[row][col] for col in col_vals]]))


def main(args):
  kv = read_input(args)
  print_matrix(kv, args)


if __name__ == "__main__":
  argparse = argparse.ArgumentParser()
  argparse.add_argument('--input', '-i', required=True)
  args = argparse.parse_args()

  main(args)
