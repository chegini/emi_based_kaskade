#!/usr/bin/env python3
"""Plot true residuals and CG step lengths from --bddcConvergenceReport CSV files."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Sample:
    time_step: int
    iteration: int
    step_length: float | None
    residual_norm: float
    rhs_norm: float
    relative_residual: float
    target: float


def read_history(path: Path) -> list[Sample]:
    required = {
        "time_step", "iteration", "step_length", "residual_norm",
        "rhs_norm", "relative_residual", "target",
    }
    samples = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing convergence CSV columns")
        for row in reader:
            try:
                sample = Sample(
                    int(row["time_step"]), int(row["iteration"]),
                    float(row["step_length"]) if row["step_length"] else None,
                    float(row["residual_norm"]), float(row["rhs_norm"]),
                    float(row["relative_residual"]), float(row["target"]),
                )
            except (TypeError, ValueError) as error:
                raise ValueError(f"{path}: malformed convergence CSV row {reader.line_num}") from error
            values = (sample.residual_norm, sample.rhs_norm, sample.relative_residual)
            if any(not math.isfinite(value) or value < 0 for value in values):
                raise ValueError(f"{path}: invalid residual/RHS norm at iteration {sample.iteration}")
            if not math.isfinite(sample.target) or sample.target <= 0:
                raise ValueError(f"{path}: invalid stopping target")
            if sample.time_step < 1 or sample.iteration != len(samples):
                raise ValueError(f"{path}: expected accepted-update indices starting at zero")
            if sample.iteration == 0 and sample.step_length is not None:
                raise ValueError(f"{path}: initial row must have no step length")
            if sample.iteration > 0 and (
                sample.step_length is None or not math.isfinite(sample.step_length)
                or sample.step_length <= 0
            ):
                raise ValueError(f"{path}: invalid CG step length")
            relative = sample.residual_norm / sample.rhs_norm if sample.rhs_norm else sample.residual_norm
            if not math.isclose(relative, sample.relative_residual, rel_tol=1e-12, abs_tol=0):
                raise ValueError(f"{path}: inconsistent relative residual")
            if samples and (sample.time_step, sample.rhs_norm, sample.target) != (
                samples[0].time_step, samples[0].rhs_norm, samples[0].target
            ):
                raise ValueError(f"{path}: timestep, RHS norm or target changes within one history")
            samples.append(sample)
    if not samples:
        raise ValueError(f"{path}: no convergence samples")
    return samples


def plot_history(path: Path, samples: list[Sample], absolute: bool, output: Path) -> None:
    # Keep CSV reading and --help usable in the solver image without Matplotlib.
    try:
        import matplotlib
    except ImportError as error:
        raise RuntimeError(
            "Plotting requires Matplotlib; install it in your plotting environment "
            "with: python3 -m pip install matplotlib"
        ) from error
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    fig, (residual_ax, alpha_ax) = plt.subplots(
        2, 1, sharex=True, figsize=(8, 6), constrained_layout=True,
        gridspec_kw={"height_ratios": [2, 1]},
    )
    first = samples[0]
    values = [sample.residual_norm if absolute else sample.relative_residual for sample in samples]
    target = first.target * first.rhs_norm if absolute and first.rhs_norm else first.target
    positive = [(sample.iteration, value) for sample, value in zip(samples, values) if value > 0]
    if positive:
        residual_ax.semilogy(
            [item[0] for item in positive], [item[1] for item in positive],
            "o-", color="tab:blue", markersize=3, label="True global residual",
        )
    else:
        residual_ax.set_yscale("log")
        residual_ax.set_ylim(target / 10, target * 10)
    residual_ax.axhline(target, color="tab:red", linestyle="--", label="Stopping target")
    if absolute and first.rhs_norm:
        residual_ax.axhline(
            first.rhs_norm, color="0.45", linestyle=":", label=r"Fixed $\|b\|_2$",
        )
    zero_at = [str(sample.iteration) for sample, value in zip(samples, values) if value == 0]
    if zero_at:
        residual_ax.text(
            0.02, 0.04, "Zero residual at update(s) " + ", ".join(zero_at)
            + " (omitted from log axis)", transform=residual_ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9},
        )
    if absolute or not first.rhs_norm:
        residual_ax.set_ylabel(r"Absolute residual $\|b-A\Delta u\|_2$")
    else:
        residual_ax.set_ylabel(r"Relative residual $\|b-A\Delta u\|_2 / \|b\|_2$")
    residual_ax.set_title(f"BDDC convergence - time step {first.time_step}\n{path.parent.name}")
    residual_ax.legend(fontsize=9)

    updates = [sample for sample in samples if sample.step_length is not None]
    if updates:
        alpha_ax.plot(
            [sample.iteration for sample in updates], [sample.step_length for sample in updates],
            "o-", color="tab:orange", markersize=3,
        )
    else:
        alpha_ax.text(0.5, 0.5, "No accepted CG updates", ha="center", va="center", transform=alpha_ax.transAxes)
    alpha_ax.set_ylabel(r"Step length $\alpha$")
    alpha_ax.set_xlabel("Accepted CG updates (0 = initial correction after interior elimination)")
    alpha_ax.set_xlim(-0.25, max(1, samples[-1].iteration) + 0.25)
    alpha_ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for axis in (residual_ax, alpha_ax):
        axis.grid(True, which="both", alpha=0.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_files", type=Path, nargs="+")
    parser.add_argument("--absolute", action="store_true", help="plot absolute residual and fixed RHS norm")
    parser.add_argument("--format", choices=("png", "pdf", "svg"), default="png")
    parser.add_argument("--output-dir", type=Path, help="default: beside each CSV")
    args = parser.parse_args()
    try:
        for path in args.csv_files:
            samples = read_history(path)
            suffix = "-absolute" if args.absolute else ""
            output = (args.output_dir or path.parent) / f"{path.stem}{suffix}.{args.format}"
            plot_history(path, samples, args.absolute, output)
            print(output)
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
