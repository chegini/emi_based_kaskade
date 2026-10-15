#!/usr/bin/env python3
"""Compare raw, zigzag, and bit-length histograms for all lifting modes."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


CASES = [
    ("No lifting", "no-graph-lifting", "#c45a48"),
    ("Lifting after quantization", "graph-lifting-after-quantization", "#3568a8"),
]


def find_case(root: Path, prefix: str) -> Path:
    for directory in sorted(root.glob(prefix + "_bits-*_histogram"), reverse=True):
        direct = directory / "restrict_raw.csv"
        nested = directory / "histograms" / "restrict_raw.csv"
        if direct.exists():
            return directory
        if nested.exists():
            return directory / "histograms"
    raise FileNotFoundError(f"no histogram directory for {prefix} below {root}")


def read(path: Path) -> list[tuple[int, float]]:
    with path.open(newline="") as stream:
        rows = []
        for row in csv.DictReader(stream):
            if not row.get("symbol", "").strip().isdigit():
                continue
            rows.append((int(row["symbol"]), float(row["probability"])))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--top", type=int, default=30,
                        help="number of most probable symbols for raw/zigzag plots")
    args = parser.parse_args()

    for stage in ("raw", "zigzag", "bitlength"):
        data: dict[tuple[str, str], list[tuple[int, float]]] = {}
        for label, prefix, _ in CASES:
            directory = find_case(args.root, prefix)
            for phase in ("restrict", "prolongate"):
                values = read(directory / f"{phase}_{stage}.csv")
                if stage != "bitlength" and args.top > 0:
                    values = sorted(values, key=lambda item: item[1], reverse=True)[:args.top]
                    values.sort(key=lambda item: item[0])
                data[(label, phase)] = values

        fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=False)
        for axis, phase in zip(axes, ("restrict", "prolongate")):
            symbols = sorted({symbol for label, _, _ in CASES
                              for symbol, _ in data[(label, phase)]})
            positions = {symbol: index for index, symbol in enumerate(symbols)}
            offsets = [-0.18, 0.18]
            for index, (label, _, color) in enumerate(CASES):
                values = dict(data[(label, phase)])
                axis.bar([positions[symbol] + offsets[index] for symbol in symbols],
                         [values.get(symbol, 0.0) for symbol in symbols],
                         width=0.24, color=color, alpha=0.86, label=label)
            axis.set_title("Restriction" if phase == "restrict" else "Prolongation")
            axis.set_xlabel("Symbol")
            axis.set_ylabel("Probability")
            axis.grid(axis="y", alpha=0.25)
            axis.set_axisbelow(True)
            if len(symbols) > 40:
                step = max(1, len(symbols) // 20)
                axis.set_xticks(range(0, len(symbols), step))
                axis.set_xticklabels([str(symbols[i]) for i in range(0, len(symbols), step)],
                                     rotation=45, ha="right")
            else:
                axis.set_xticks(range(len(symbols)))
                axis.set_xticklabels([str(symbol) for symbol in symbols], rotation=45,
                                     ha="right")
        axes[1].legend(frameon=False, loc="upper right")
        suffix = "top" if stage != "bitlength" else "all"
        fig.suptitle(f"{args.title}: {stage} symbols ({suffix})")
        fig.tight_layout()
        output = args.output_prefix.with_name(
            f"{args.output_prefix.name}-{stage}")
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output.with_suffix(".png"), dpi=220)
        fig.savefig(output.with_suffix(".pdf"))
        plt.close(fig)


if __name__ == "__main__":
    main()
