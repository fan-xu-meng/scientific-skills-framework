#!/usr/bin/env python3
"""Check the approximate independent-point budget for an EXAFS fit window."""

from __future__ import annotations

import argparse
import math


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kmin", type=float, required=True)
    parser.add_argument("--kmax", type=float, required=True)
    parser.add_argument("--rmin", type=float, required=True)
    parser.add_argument("--rmax", type=float, required=True)
    parser.add_argument("--free-parameters", type=int, required=True)
    args = parser.parse_args()
    values = (args.kmin, args.kmax, args.rmin, args.rmax)
    if not all(math.isfinite(v) for v in values) or args.kmax <= args.kmin or args.rmax <= args.rmin or args.free_parameters < 0:
        parser.error("use finite increasing k/R ranges and a nonnegative parameter count")
    nind = 2 * (args.kmax - args.kmin) * (args.rmax - args.rmin) / math.pi
    print(f"Approximate Nind: {nind:.3f}")
    print(f"Free parameters: {args.free_parameters}")
    print(f"Budget status: {'overparameterized' if args.free_parameters >= nind else 'within approximate budget'}")
    print("This is a necessary complexity check, not evidence that the structural model is correct.")


if __name__ == "__main__":
    main()
