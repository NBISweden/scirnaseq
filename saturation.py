#!/usr/bin/env python3
"""
Read in a BAM file and create necessary data for a saturation plot.

- CR and UR tags are used because these can be output by STAR solo even without
  needing to sort the output BAM. Since we use simulated cell barcodes, this
  should mostly be ok - UMIs are not simulated though, so this introduces an
  error.

- The (UR, CR, GX) triple determines uniqueness.
"""
import sys
from argparse import ArgumentParser
from pathlib import Path
import subprocess

from pysam import AlignmentFile


def main():
    parser = ArgumentParser()
    parser.add_argument("--barcodes", type=Path, help="File with allowed cell barcodes")
    parser.add_argument("bam")
    args = parser.parse_args()


    if args.barcodes:
        allowed_cell_barcodes = set(read_barcodes(args.barcodes))
        n_cells = len(allowed_cell_barcodes)
        print("Read file with", n_cells, "allowed cell barcodes", file=sys.stderr)
    else:
        allowed_cell_barcodes = None
        n_cells = 0

    combinations = set()
    umis = set()
    n_not_allowed = 0

    print("n", "n_per_cell", "saturation", sep="\t")
    print(0, 0, 0, sep="\t")
    n = 0
    for record in record_iterator(args.bam):
        if record.is_secondary or record.mapping_quality < 255:
            continue
        cell_barcode = record.get_tag("CR")
        if allowed_cell_barcodes is not None and cell_barcode not in allowed_cell_barcodes:
            n_not_allowed += 1
            continue
        umi = record.get_tag("UR")
        gx = record.get_tag("GX")
        if gx == "-":
            continue

        n += 1
        n_per_cell = n / n_cells if n_cells is not None else 0
        combinations.add((cell_barcode, umi))
        if n % 10000 == 0:
            saturation = 1 - (len(combinations) / n)
            print(f"{saturation:6.4f} {int(saturation*100)*'*'}", file=sys.stderr)
            print(n, f"{saturation:.4f}", sep="\t")

    saturation = 1 - (len(combinations) / n)
    print(f"{saturation:6.4f} {int(saturation*100)*'*'}", file=sys.stderr)
    print(n, n_per_cell, f"{saturation:.4f}", sep="\t")
    print(file=sys.stderr)
    print(n_not_allowed, "records with disallowed cell barcodes skipped", file=sys.stderr)
    print(n, "remaining records", file=sys.stderr)
    print(len(combinations), "unique cell barcode/UMI combinations", file=sys.stderr)
    print("saturation:", 1 - (len(combinations) / n), file=sys.stderr)


def read_barcodes(path):
    with open(path) as f:
        for line in f:
            yield line.strip()


def record_iterator(path):
    with subprocess.Popen(["samtools", "collate", "-u", "-O", path], stdout=subprocess.PIPE) as samtools:
        with AlignmentFile(samtools.stdout) as af:
            yield from af


if __name__ == "__main__":
    main()
