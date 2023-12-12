#!/usr/bin/env python3
"""
Simulate a cell barcode
"""
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from itertools import product
import sys
import re

import dnaio
from tinyalign import hamming_distance


@dataclass
class RecordHeader:
    i1: str
    ligation_index: str
    rt_index: str
    umi: str


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--r1", "-o", help="R1 output (simulated cell barcode)", type=Path
    )
    parser.add_argument("--r2", "-p", help="R2 output (cDNA)", type=Path)
    parser.add_argument(
        "--list",
        metavar="FILE",
        help="Write all possible emulated barcodes to FILE",
        type=Path,
    )
    parser.add_argument(
        "--p7-stats",
        metavar="FILE",
        help="Write P7 index (I1) statistics to FILE",
    )
    parser.add_argument("ligation_indices_fasta", type=Path)
    parser.add_argument("rt_indices_fasta", type=Path)
    parser.add_argument("p7_indices_fasta", type=Path)
    parser.add_argument("fastq", nargs="+", type=Path)
    args = parser.parse_args()

    # All ligation indices end with "T".
    # We append an "A" to those that only have 9 nucleotides.
    ligation_indices = make_all_same_length(read_fasta(args.ligation_indices_fasta))
    rt_indices = read_fasta(args.rt_indices_fasta)
    p7_indices = read_fasta(args.p7_indices_fasta)

    if args.list:
        write_list(args.list, rt_indices, ligation_indices, p7_indices)
        print("Wrote", args.list, file=sys.stderr)

    p7_index_mismatches = []

    r1_path = args.r1
    r2_path = args.r2
    with dnaio.open(r1_path, r2_path, mode="w") as outf:
        for fastq_path in args.fastq:
            #sample = Path(fastq_path.stem).stem  # S1, ..., S96
            sample = re.search("_(S[0-9]+)_", fastq_path.stem)[1]
            lane = re.search("_(L[0-9]+)_", fastq_path.stem)[1]
            p7_index = p7_indices[sample]
            mismatches = [0] * len(p7_index)
            with dnaio.open(fastq_path) as inf:
                for record in inf:
                    header = parse(record.name)
                    mismatches[hamming_distance(p7_index, header.i1)] += 1
                    record.name = record.id
                    barcode_record = record[:]
                    components = (
                        rt_indices[header.rt_index],
                        ligation_indices[header.ligation_index],
                        p7_index,
                        header.umi,
                    )
                    barcode_record.sequence = "".join(components)
                    barcode_record.qualities = "F" * len(barcode_record.sequence)
                    outf.write(barcode_record, record)

            p7_index_mismatches.append([sample, lane, p7_index] + mismatches)


    if args.p7_stats:
        write_mismatch_stats(args.p7_stats, p7_index_mismatches)
        print("Wrote", args.p7_stats, file=sys.stderr)


def read_fasta(path):
    """
    Read a FASTA file and return a dict that maps a record name to a sequence
    """
    with dnaio.open(path) as f:
        return {record.id: record.sequence.upper() for record in f}


def make_all_same_length(d):
    target_length = max(len(seq) for seq in d.values())

    return {name: seq.ljust(target_length, "A") for name, seq in d.items()}


def parse(header):
    fields = header.split(" ")
    i1 = fields[1].split(":")[3]
    assert fields[2].startswith("ligation_index=")
    ligation_index = fields[2].split("=")[1]
    assert fields[3].startswith("umi=")
    umi = fields[3].split("=")[1]
    assert fields[4].startswith("rt_index=")
    rt_index = fields[4].split("=")[1]

    return RecordHeader(i1=i1, ligation_index=ligation_index, rt_index=rt_index, umi=umi)


def write_list(path, *indices_dicts):
    with open(path, "w") as f:
        for index_keys in product(*indices_dicts):
            # name = "_".join(index_keys)
            seq = "".join(
                index_dict[name] for index_dict, name in zip(indices_dicts, index_keys)
            )
            print(seq, file=f)


def write_mismatch_stats(path, stats):
    with open(path, "w") as f:
        print("sample", "lane", "sequence", "total", "perfect", "one_mismatch", "two_mismatches", sep="\t", file=f)
        for row in stats:
            print(*row[0:3], sum(row[3:]), *row[3:6], sep="\t", file=f)


if __name__ == "__main__":
    main()
