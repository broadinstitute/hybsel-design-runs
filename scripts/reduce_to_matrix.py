#!/bin/python3

import argparse
import itertools
from collections import defaultdict


def read_input(args):
  kv = defaultdict(dict)
  with open(args.input) as f:
    for line in f:
      ls = line.rstrip().split()
      if len(ls) < 3:
        continue
      x, y, val = \
          float(ls[args.row_index]), float(ls[args.col_index]), float(ls[args.val_index])
      if int(x) == x:
        x = int(x)
      if int(y) == y:
        y = int(y)
      if int(val) == val:
        val = int(val)
      if args.require_index and args.require_val:
        if float(ls[args.require_index]) != args.require_val:
          # skip this entry
          continue
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
  argparse.add_argument('--row_index', '-r', type=int, default=0)
  argparse.add_argument('--col_index', '-c', type=int, default=1)
  argparse.add_argument('--val_index', '-v', type=int, default=2)
  argparse.add_argument('--require_index', type=int)
  argparse.add_argument('--require_val', type=float)
  args = argparse.parse_args()

  main(args)
