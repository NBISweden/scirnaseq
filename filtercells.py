#!/usr/bin/env python3
"""
Simple cell filtering (keeps only cells with a certain minimum UMI count)
"""
from argparse import ArgumentParser
from pathlib import Path
import sys
import shutil

import scipy
from xopen import xopen


def main():
    parser = ArgumentParser()
    parser.add_argument("--umis", type=int, default=200, help="Minimum UMI count")
    parser.add_argument("raw", type=Path, help="Input directory with matrix.mtx, features.tsv and barcodes.tsv file")
    parser.add_argument("target", type=Path, help="Output directory (must not exist)")
    args = parser.parse_args()

    if args.target.exists():
        parser.error(f"Target {args.target} already exists")

    log("Reading barcodes.tsv(.gz) ...")
    try:
        with xopen(args.raw / "barcodes.tsv") as f:
            barcodes = list(f)
        barcodes_name = "barcodes.tsv"
    except FileNotFoundError:
        with xopen(args.raw / "barcodes.tsv.gz") as f:
            barcodes = list(f)
        barcodes_name = "barcodes.tsv.gz"

    log("Reading matrix.mtx(.gz) ...")
    try:
        matrix = scipy.io.mmread(args.raw / "matrix.mtx.gz")
        matrix_name = "matrix.mtx.gz"
    except FileNotFoundError:
        matrix = scipy.io.mmread(args.raw / "matrix.mtx")
        matrix_name = "matrix.mtx"

    log(f"Input matrix: {matrix.shape[1]} cell barcodes; {matrix.shape[0]} features")

    log("Filtering ...")
    matrix = matrix.tocsr()
    umi_counts = matrix.sum(axis=0).A1
    mask = umi_counts >= args.umis

    filtered_matrix = matrix[:, mask]
    filtered_barcodes = [barcode for barcode, use_it in zip(barcodes, mask) if use_it]

    log(f"Filtered matrix: {filtered_matrix.shape[1]} cell barcodes; {filtered_matrix.shape[0]} features")

    log("Writing output ...")
    args.target.mkdir()
    with xopen(args.target / barcodes_name, mode="w") as f:
        for barcode in filtered_barcodes:
            f.write(barcode)

    try:
        shutil.copy(args.raw / "features.tsv.gz", args.target / "features.tsv.gz")
    except FileNotFoundError:
        shutil.copy(args.raw / "features.tsv", args.target / "features.tsv")

    scipy.io.mmwrite(args.target / matrix_name, filtered_matrix)


def log(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


if __name__ == "__main__":
    main()
