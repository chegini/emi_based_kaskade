#!/usr/bin/env python3
"""Collect and validate solver logs for container checks, not scaling timings."""
import argparse
import json
import math
from pathlib import Path
import re
import sys

NUMBER = r"[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?"


def collect(text, expected_steps, expected_ranks, exit_code, expected_tolerance=1e-7):
    if not math.isfinite(expected_tolerance) or expected_tolerance <= 0:
        raise ValueError("Expected tolerance must be finite and positive")
    errors = []
    rows = re.findall(
        rf"^BDDC global residual step (\d+): abs=({NUMBER}) rel=({NUMBER}) "
        rf"target=({NUMBER}) status=(\S+) iterations=(\d+)\s*$", text, re.MULTILINE)
    if exit_code:
        errors.append(f"MPI/application exit code {exit_code}")
    if not re.search(r"^EMI model completed\s*$", text, re.MULTILINE):
        errors.append("Missing normal completion marker")
    if [int(row[0]) for row in rows] != list(range(1, expected_steps+1)):
        errors.append("Missing, duplicate or out-of-order true-residual summaries")
    for step, absolute, relative, target, status, iterations in rows:
        values = [float(absolute), float(relative), float(target)]
        if not all(math.isfinite(value) and value >= 0 for value in values):
            errors.append(f"Nonfinite/negative residual in step {step}")
        if status != "converged" or float(target) != expected_tolerance or float(relative) > expected_tolerance:
            errors.append(f"Step {step} did not meet the true relative tolerance {expected_tolerance:g}")

    def single(pattern, label, cast=float):
        matches = re.findall(pattern, text, re.MULTILINE)
        if len(matches) != 1:
            errors.append(f"Expected one {label}; found {len(matches)}")
            return None
        value = cast(matches[0])
        if not math.isfinite(value):
            errors.append(f"Nonfinite {label}")
            return None
        return value

    cells = single(r"^cells:\s+(\d+)\s*$", "cell count", int)
    dofs = single(r"^dofs:\s+(\d+)\s*$", "DOF count", int)
    subdomains = single(r"^BDDC subdomains:\s+(\d+),.*$", "subdomain count", int)
    shared = single(r"^BDDC subdomains:\s+\d+, shared dofs:\s+(\d+)\s*$",
                    "partition shared DOF count", int)
    coarse = re.findall(r"^BDDC coarse interfaces: corners=(\d+), edges=(\d+), faces=(\d+)\s*$",
                        text, re.MULTILINE)
    if len(coarse) != 1:
        errors.append("Expected one set of coarse-interface counts")
    ranks = single(r"^Profile \(wall seconds, EMI end-to-end phases; (\d+) rank\(s\)\):$",
                   "phase-profile rank count", int)
    if ranks != expected_ranks:
        errors.append("Actual rank count differs from requested count")
    raw = single(r"^\s*raw sent bytes:\s+(\d+)\s*$", "raw sent bytes", int)
    wire = single(r"^\s*transmitted sent bytes:\s+(\d+)\s*$", "transmitted sent bytes", int)
    wall = single(rf"^total cpu-time:\s+({NUMBER})s wall,.*$", "application wall time")
    wrapper = single(rf"^\s*solver wrapper \(including setup\)\s+{NUMBER}\s+{NUMBER}\s+({NUMBER})\s*$",
                     "solver-wrapper maximum")
    return {
        "status": "failed" if errors else "converged",
        "errors": errors,
        "validation_only": True,
        "exit_code": exit_code,
        "requested_relative_tolerance": expected_tolerance,
        "cells": cells,
        "dofs": dofs,
        "subdomains": subdomains,
        "partition_shared_dofs": shared,
        "coarse_interfaces": dict(zip(("corners","edges","faces"),map(int,coarse[0])))
                             if len(coarse)==1 else None,
        "mpi_ranks": ranks,
        "iterations_per_step": [int(row[5]) for row in rows],
        "relative_residuals": [float(row[2]) if math.isfinite(float(row[2])) else None
                               for row in rows],
        "raw_sent_bytes": raw,
        "transmitted_sent_bytes": wire,
        "reported_transfer_savings_percent": 100*(1-wire/raw) if raw and wire is not None else None,
        "wall_seconds": wall,
        "solver_wrapper_max_seconds": wrapper,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--exit-code", type=int, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-7)
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error("--tolerance must be finite and positive")
    result = collect(args.log.read_text(errors="replace"), args.steps, args.ranks, args.exit_code,
                     args.tolerance)
    output = args.log.with_name("metrics.json")
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    if result["errors"]:
        print("FAILED: "+"; ".join(result["errors"]), file=sys.stderr)
        return 1
    savings=result["reported_transfer_savings_percent"]
    savings_text=f"{savings:.2f}%" if savings is not None else "n/a (zero raw bytes)"
    print(f"Converged: iterations={result['iterations_per_step']}, "
          f"max_rel={max(result['relative_residuals']):.6g}, "
          f"reported_transfer_savings={savings_text}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
