#!/usr/bin/env python3
"""
Write out a list of all possibel (simulated) cell barcodes
"""
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from itertools import product
import sys
import re

import dnaio


def main():
    parser = ArgumentParser()
    parser.add_argument("ligation_indices_fasta", type=Path)
    parser.add_argument("rt_indices_fasta", type=Path)
    parser.add_argument("p7_indices_fasta", type=Path)
    args = parser.parse_args()

    # All ligation indices end with "T".
    # We append an "A" to those that only have 9 nucleotides.
    ligation_indices = make_all_same_length(read_fasta(args.ligation_indices_fasta))
    rt_indices = read_fasta(args.rt_indices_fasta)
    p7_indices = read_fasta(args.p7_indices_fasta)

    indices_dicts = rt_indices, ligation_indices, p7_indices
    n = 0
    for index_keys in product(*indices_dicts):
        # name = "_".join(index_keys)
        seq = "".join(
            index_dict[name] for index_dict, name in zip(indices_dicts, index_keys)
        )
        print(seq)
        n += 1
    print("Wrote", n, "possible barcodes", file=sys.stderr)


def read_fasta(path):
    """
    Read a FASTA file and return a dict that maps a record name to a sequence
    """
    with dnaio.open(path) as f:
        return {record.id: record.sequence.upper() for record in f}


def make_all_same_length(d):
    target_length = max(len(seq) for seq in d.values())

    return {name: seq.ljust(target_length, "A") for name, seq in d.items()}


if __name__ == "__main__":
    main()
