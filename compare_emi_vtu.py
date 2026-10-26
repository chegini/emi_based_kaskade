#!/usr/bin/env python3
"""Compare final EMI solutions stored in ASCII VTK .vtu files."""

from __future__ import annotations

import argparse
import math
import xml.etree.ElementTree as ET
from pathlib import Path


def read_array(path: Path, name: str) -> list[float]:
    root = ET.parse(path).getroot()
    for array in root.iter("DataArray"):
        if array.attrib.get("Name") == name:
            text = array.text or ""
            return [float(value) for value in text.split()]
    raise ValueError(f"DataArray {name!r} not found in {path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare the scalar solution array u in two ASCII VTU files."
    )
    parser.add_argument("compressed", type=Path, help="compressed solution VTU")
    parser.add_argument("reference", type=Path, help="uncompressed reference VTU")
    args = parser.parse_args()

    compressed_coordinates = read_array(args.compressed, "Coordinates")
    reference_coordinates = read_array(args.reference, "Coordinates")
    compressed = read_array(args.compressed, "u")
    reference = read_array(args.reference, "u")

    if len(compressed_coordinates) != len(reference_coordinates):
        raise ValueError("coordinate arrays have different lengths")
    if any(a != b for a, b in zip(compressed_coordinates, reference_coordinates)):
        raise ValueError("VTU files do not have identical coordinates/order")
    if len(compressed) != len(reference):
        raise ValueError("solution arrays have different lengths")

    errors = [a - b for a, b in zip(compressed, reference)]
    abs_errors = [abs(value) for value in errors]
    reference_abs = [abs(value) for value in reference]
    l2_error = math.sqrt(sum(value * value for value in errors))
    l2_reference = math.sqrt(sum(value * value for value in reference))
    max_error = max(abs_errors, default=0.0)
    max_reference = max(reference_abs, default=0.0)
    mean_error = sum(abs_errors) / len(abs_errors) if abs_errors else 0.0

    print(f"compressed:          {args.compressed}")
    print(f"reference:           {args.reference}")
    print(f"solution values:     {len(compressed)}")
    print("same coordinates:    yes")
    print(f"max absolute error:  {max_error:.16e}")
    print(f"L2 error:            {l2_error:.16e}")
    print(f"L2 reference:        {l2_reference:.16e}")
    print(f"relative L2 error:   {l2_error / l2_reference:.16e}")
    print(f"relative max error:  {max_error / max_reference:.16e}")
    print(f"mean absolute error: {mean_error:.16e}")


if __name__ == "__main__":
    main()
