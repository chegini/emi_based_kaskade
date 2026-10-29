#!/usr/bin/env python3
"""Compare normalized old/scratch BDDC local-space and coarse-space dumps."""

import csv
import math
import sys
from pathlib import Path


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def keyed_values(path, key_fields, value_field):
    result = {}
    for row in rows(path):
        key = tuple(row[field] for field in key_fields)
        result[key] = result.get(key, 0.0) + float(row[value_field])
    return result


def read_global_matrix(path):
    result = {}
    for row in rows(path):
        key = (row["row"], row["col"])
        result[key] = result.get(key, 0.0) + float(row["value"])
    return result


def compare_sets(label, a, b):
    only_a = sorted(a - b)
    only_b = sorted(b - a)
    print(f"{label}: old={len(a)} scratch={len(b)} common={len(a & b)}")
    print(f"  only old={len(only_a)} only scratch={len(only_b)}")
    if only_a:
        print(f"  first only old: {only_a[:5]}")
    if only_b:
        print(f"  first only scratch: {only_b[:5]}")


def compare_values(label, a, b):
    ka, kb = set(a), set(b)
    compare_sets(label + " pattern", ka, kb)
    errors = [a[k] - b[k] for k in ka & kb]
    max_error = max((abs(value) for value in errors), default=0.0)
    l2_error = math.sqrt(sum(value * value for value in errors))
    print(f"  common-entry max abs difference={max_error:.17g}; L2={l2_error:.17g}")
    if max_error:
        worst = max(ka & kb, key=lambda key: abs(a[key] - b[key]))
        print(f"  largest common-entry mismatch at {worst}: old={a[worst]:.17g} scratch={b[worst]:.17g}")


def main():
    if len(sys.argv) not in (3, 6) or (len(sys.argv) == 6 and sys.argv[3] != "--global"):
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} OLD_BDDC_LOCAL_DATA_DIR SCRATCH_BDDC_LOCAL_DATA_DIR [--global OLD_LHS.csv SCRATCH_LHS.csv]")
    old, scratch = (Path(arg) for arg in sys.argv[1:3])

    old_dofs = keyed_values(old / "local_dofs_weights.csv", ["tag", "global_dof"], "weight")
    new_dofs = keyed_values(scratch / "local_dofs_weights.csv", ["tag", "global_dof"], "weight")
    compare_values("local memberships and weights", old_dofs, new_dofs)

    old_order = {(r["tag"], r["local_dof"], r["global_dof"])
                 for r in rows(old / "local_dofs_weights.csv")}
    new_order = {(r["tag"], r["local_dof"], r["global_dof"])
                 for r in rows(scratch / "local_dofs_weights.csv")}
    compare_sets("local-index to global-index maps", old_order, new_order)

    old_matrices = keyed_values(old / "local_matrices.csv",
                                ["tag", "row_global", "col_global"], "value")
    new_matrices = keyed_values(scratch / "local_matrices.csv",
                                ["tag", "row_global", "col_global"], "value")
    compare_values("local matrix entries", old_matrices, new_matrices)
    old_local_rhs = keyed_values(old / "local_rhs.csv", ["tag", "global_dof"], "value")
    new_local_rhs = keyed_values(scratch / "local_rhs.csv", ["tag", "global_dof"], "value")
    compare_values("local RHS by material tag and global DOF", old_local_rhs, new_local_rhs)
    if len(sys.argv) == 6:
        old_lhs = read_global_matrix(Path(sys.argv[4]))
        scratch_lhs = read_global_matrix(Path(sys.argv[5]))
        old_sum = {}
        for key, value in old_matrices.items():
            global_key = (key[1], key[2])
            old_sum[global_key] = old_sum.get(global_key, 0.0) + value
        new_sum = {}
        for key, value in new_matrices.items():
            global_key = (key[1], key[2])
            new_sum[global_key] = new_sum.get(global_key, 0.0) + value
        compare_values("old embedded local-matrix sum vs global LHS", old_lhs, old_sum)
        compare_values("scratch embedded local-matrix sum vs global LHS", scratch_lhs, new_sum)
        old_rhs = keyed_values(Path(sys.argv[4]).with_name("global_rhs.csv"), ["dof"], "value")
        scratch_rhs = keyed_values(Path(sys.argv[5]).with_name("global_rhs.csv"), ["dof"], "value")
        old_rhs_sum = keyed_values(old / "local_rhs.csv", ["global_dof"], "value")
        new_rhs_sum = keyed_values(scratch / "local_rhs.csv", ["global_dof"], "value")
        compare_values("old local-RHS sum vs global RHS", old_rhs, old_rhs_sum)
        compare_values("scratch local-RHS sum vs global RHS", scratch_rhs, new_rhs_sum)

    old_shared = keyed_values(old / "shared_dofs.csv", ["global_dof", "tag"], "weight")
    new_shared = keyed_values(scratch / "shared_dofs.csv", ["global_dof", "tag"], "weight")
    compare_values("shared-DOF memberships and weights", old_shared, new_shared)

    old_coarse = {(r["type"], r["tag"], r["global_dof"]) for r in rows(old / "coarse_constraints.csv")}
    new_coarse = {(r["type"], r["tag"], r["global_dof"]) for r in rows(scratch / "coarse_constraints.csv")}
    compare_sets("coarse constraint memberships", old_coarse, new_coarse)

    old_summary = (old / "coarse_summary.txt").read_text().splitlines()
    new_summary = (scratch / "coarse_summary.txt").read_text().splitlines()
    print(f"coarse summaries identical: {old_summary == new_summary}")
    if old_summary != new_summary:
        print(f"  old: {old_summary}")
        print(f"  scratch: {new_summary}")


if __name__ == "__main__":
    main()
