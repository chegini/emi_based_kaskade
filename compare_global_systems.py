#!/usr/bin/env python3
"""Compare global_lhs.csv/global_rhs.csv dumps from the old and scratch codes."""

import csv
import math
import sys
from pathlib import Path


def read_matrix(path):
    values = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            key = (int(row["row"]), int(row["col"]))
            values[key] = values.get(key, 0.0) + float(row["value"])
    return values


def read_vector(path):
    values = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            values[int(row["dof"])] = float(row["value"])
    return values


def compare(label, a, b):
    keys = set(a) | set(b)
    missing_a = sum(k not in a for k in keys)
    missing_b = sum(k not in b for k in keys)
    errors = [a.get(k, 0.0) - b.get(k, 0.0) for k in keys]
    max_error = max((abs(x) for x in errors), default=0.0)
    l2_error = math.sqrt(sum(x * x for x in errors))
    print(f"{label}:")
    print(f"  entries: old={len(a)} scratch={len(b)} union={len(keys)}")
    print(f"  missing: old={missing_a} scratch={missing_b}")
    print(f"  max absolute difference: {max_error:.17g}")
    print(f"  L2 difference: {l2_error:.17g}")
    if max_error:
        worst = max(keys, key=lambda k: abs(a.get(k, 0.0) - b.get(k, 0.0)))
        print(f"  largest mismatch at {worst}: old={a.get(worst, 0.0):.17g} scratch={b.get(worst, 0.0):.17g}")


def main():
    if len(sys.argv) != 3:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} OLD_GLOBAL_SYSTEM_DIR SCRATCH_GLOBAL_SYSTEM_DIR")
    old_dir, scratch_dir = map(Path, sys.argv[1:])
    compare("LHS", read_matrix(old_dir / "global_lhs.csv"),
            read_matrix(scratch_dir / "global_lhs.csv"))
    compare("RHS", read_vector(old_dir / "global_rhs.csv"),
            read_vector(scratch_dir / "global_rhs.csv"))


if __name__ == "__main__":
    main()
