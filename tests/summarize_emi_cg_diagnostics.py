#!/usr/bin/env python3
"""Summarize read-only EMI CG observations; this does not certify convergence."""
import argparse
import csv
import math
from pathlib import Path


def finite(row, key):
    try:
        value = float(row[key])
        return value if math.isfinite(value) else None
    except (KeyError, ValueError):
        return None


def summarize(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    restriction = [r for r in rows if r['stage'] == 'restriction']
    search = [r for r in rows if r['stage'] == 'search']
    def maximum(samples, key):
        values = [v for r in samples if (v := finite(r,key)) is not None]
        return f'{max(values):.4e}' if values else 'undefined'
    print(path)
    for p in sorted({int(r['restriction_pass']) for r in restriction}):
        selected = [r for r in restriction if int(r['restriction_pass']) == p]
        print(f"  restriction pass {p}: max relative error = "
              + maximum(selected,'restriction_error_rel'))
    print("  max relative rho discrepancy = "+maximum(search,'rho_error_rel'))
    print("  max adjacent-direction conjugacy defect = "+maximum(search,'conjugacy_previous'))
    print("  max beta = "+maximum(search,'beta'))
    print("  max candidate alpha = "+maximum(search,'alpha_candidate'))
    print("  rejected search observations = "+str(sum(r['applied']=='0' for r in search)))
    for key in ['rho_internal','rho_true','curvature_internal','curvature_global']:
        count = sum(v <= 0 for r in search if (v := finite(r,key)) is not None)
        print(f"  nonpositive {key}: {count}")
    print("  near-zero/invalid reference rho: "
          +str(sum(r.get('rho_true_near_zero') == '1' for r in search)))
    print("  max solution change between restriction observations = "
          +maximum(rows,'solution_change_since_restriction'))
    for key in ['solution_jump','z_jump','q_jump']:
        print(f"  max {key} = "+maximum(rows,key))
    print()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path,
                        help='CSV files, individual case directories, or diagnostic output roots')
    args = parser.parse_args()
    found = set()
    for p in args.paths:
        if p.is_file():
            found.add(p)
        else:
            found.update(p.glob('bddc-cg-diagnostics-step*.csv'))
            found.update(p.glob('*/bddc-cg-diagnostics-step*.csv'))
    if not found:
        parser.error('No CG diagnostic CSVs found; run with --cg-diagnostics')
    for p in sorted(found):
        summarize(p)
    print('These are pre-update diagnostics. Use the convergence CSV for the stopping result.')
    print('Blank/undefined ratios are not zero errors. First-direction conjugacy is undefined.')
    print('Instrumentation timings include extra collectives and matrix-vector products.')


if __name__ == '__main__':
    main()
