#!/usr/bin/env python3
"""
Derive P7 indices from FASTQ headers. Output is written in FASTA format to stdout.
"""
from argparse import ArgumentParser
from collections import Counter
from pathlib import Path
from itertools import islice
import sys
import re

import dnaio


def main():
    parser = ArgumentParser()
    parser.add_argument("fastq", nargs="+", type=Path)
    args = parser.parse_args()

    p7_indices = {}
    for path in args.fastq:
        # Derive sample name from FASTQ name
        #sample = Path(fastq_path.stem).stem  # S1, ..., S96
        sample = re.search("_(S[0-9]+)_", path.stem)[1]
        #if lane_match := re.search("_(L[0-9]+)_", fastq_path.stem):
        #    lane = lane_match[1]
        #else:
        #    lane = None
        counter = Counter()
        with dnaio.open(path) as infile:
            for record in islice(infile, 10000):
                i1 = record.name.split(" ")[1].split(":")[3].upper()
                counter[i1] += 1
        mostcommon = counter.most_common(2)
        if len(mostcommon) == 2 and mostcommon[0][1] < mostcommon[1][1] * 10:
            raise ValueError(f"Cannot derive correct i7 sequence from FASTQ file {path}")
        i1_sequence = mostcommon[0][0]
        if sample in p7_indices and p7_indices[sample] != i1_sequence:
            raise ValueError(f"Found conflicting i1/P7 sequences for sample {sample}")
        p7_indices[sample] = i1_sequence

    for sample, sequence in p7_indices.items():
        print(f">{sample}", sequence, sep="\n")


if __name__ == "__main__":
    main()
