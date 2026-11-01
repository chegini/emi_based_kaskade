#!/usr/bin/env python3
"""Offline buffer ranges and explicit wire-cost models; does not run a solver."""
import argparse
import csv
import json
import math
from pathlib import Path
import re


def align8(n):
    return (n+7)//8*8


def packet_cost(indices, excluded, group_of):
    """Proposed contiguous-run descriptors: uint32 count + two float32 values.

    Receiver already knows packet DOF order. Repeated appearances of one
    producer group reuse its scale, even across packets. No per-value group IDs.
    Full packet has uint32 number-of-runs and existing 8-byte stream lengths.
    This is a proposed layout model, NOT the existing Huffman implementation.
    """
    groups = [group_of[i] for i in indices if i not in excluded]
    runs = sum(k == 0 or g != groups[k-1] for k, g in enumerate(groups))
    n = len(groups)
    exact = sum(i in excluded for i in indices)
    # Existing scale/min (8), stream lengths (8), exact doubles, alignment.
    old = align8(16+2*n) + 8*exact
    new = align8(12+12*runs+(12*n+7)//8) + 8*exact
    # Same coded stream hypothetical: only metadata changes, before padding.
    return dict(values=n, exact=exact, runs=runs,
                metadata_extra=12*runs-4, packed16=old, packed12=new)


def percentile(values, p):
    if not values:
        return None
    a = sorted(values)
    return a[min(len(a)-1, int(p*(len(a)-1)))]


def analyze(root, sizes):
    snapshots = sorted(root.glob('sub*-call*.csv'))
    if not snapshots:
        raise ValueError('No restriction captures found: '+str(root))
    groups_out, buffers_out = [], []
    totals = {g: dict(snapshots=0, values=0, at16bound=0, zero_values=0,
                     packets=0, runs=0, metadata_extra=0,
                     packed16=0, packed12=0, ratios=[]) for g in sizes}
    topology_cache = {}
    for path in snapshots:
        match = re.fullmatch(r'sub(-?\d+)-call(\d+)\.csv', path.name)
        sub, call = map(int, match.groups())
        if sub not in topology_cache:
            packets, excluded = {}, set()
            with (root/f'sub{sub}-topology.csv').open() as f:
                for row in csv.DictReader(f):
                    i, neighbor = int(row['index']), int(row['neighbor'])
                    packets.setdefault(neighbor, []).append(i)
                    if int(row['excluded']):
                        excluded.add(i)
            topology_cache[sub] = packets, excluded
        packets, excluded = topology_cache[sub]
        vals = {}
        with path.open() as f:
            for row in csv.DictReader(f):
                if int(row['component']) != 0:
                    raise ValueError('This EMI analyzer requires scalar (m=1) buffers')
                i = int(row['index'])
                if bool(int(row['excluded'])) != (i in excluded):
                    raise ValueError('Capture/topology exclusion mismatch')
                if i in vals:
                    raise ValueError('Duplicate buffer entry')
                vals[i] = float(row['value'])
        if not all(math.isfinite(x) for x in vals.values()):
            raise ValueError('Non-finite captured value')
        if set(vals) != {i for indices in packets.values() for i in indices}:
            raise ValueError('Buffer/topology membership mismatch')
        # A producer DOF gets one group, reused for every neighboring copy.
        # Compressible ordinal excludes exact-corner entries.
        active = sorted(i for i in vals if i not in excluded)
        if not active:
            continue
        lo, hi = min(vals[i] for i in active), max(vals[i] for i in active)
        span = hi-lo
        maxabs = max(abs(lo), abs(hi))
        buffers_out.append(dict(subdomain=sub, invocation=call, count=len(active),
                                minimum=lo, maximum=hi, span=span, maxabs=maxabs,
                                zeros=sum(vals[i] == 0 for i in active)))
        for size in sizes:
            t = totals[size]
            t['snapshots'] += 1
            t['values'] += len(active)
            t['zero_values'] += sum(vals[i] == 0 for i in active)
            group_of = {i: ordinal//size for ordinal, i in enumerate(active)}
            for start in range(0, len(active), size):
                ids = active[start:start+size]
                v = [vals[i] for i in ids]
                gmin, gmax = min(v), max(v)
                gspan = gmax-gmin
                gabs = max(abs(gmin), abs(gmax))
                # Ideal quantization step bound, not actual float32-wire error.
                # Compare symmetric zero-preserving 12-bit half-step to the
                # original affine 16-bit half-step on the SAME captured buffer.
                ratio = (gabs/2047)/(span/65535) if span else None
                if ratio is not None:
                    t['ratios'].extend([ratio]*len(ids))
                    if ratio <= 1:
                        t['at16bound'] += len(ids)
                groups_out.append(dict(subdomain=sub, invocation=call, group_size=size,
                    group=start//size, count=len(v), minimum=gmin, maximum=gmax,
                    span=gspan, maxabs=gabs, zeros=sum(x == 0 for x in v),
                    span_fraction=gspan/span if span else None,
                    maxabs_fraction=gabs/maxabs if maxabs else None,
                    symmetric12_step_over_affine16=ratio))
            for indices in packets.values():
                cost = packet_cost(indices, excluded, group_of)
                t['packets'] += 1
                for key in ['runs', 'metadata_extra', 'packed16', 'packed12']:
                    t[key] += cost[key]
    for t in totals.values():
        ratios = t.pop('ratios')
        t['symmetric12_step_ratio_p50'] = percentile(ratios, .5)
        t['symmetric12_step_ratio_p95'] = percentile(ratios, .95)
        t['fraction_values_at_affine16_bound'] = (
            t['at16bound']/len(ratios) if ratios else None)
        t['packed_model_saving_percent'] = (
            100*(1-t['packed12']/t['packed16']) if t['packed16'] else None)
    return buffers_out, groups_out, totals


def write_csv(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('capture', type=Path)
    p.add_argument('--groups', nargs='+', type=int, default=[32,64,128,256,512])
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    if any(g < 1 for g in args.groups):
        p.error('Group sizes must be positive')
    buffers, groups, totals = analyze(args.capture, args.groups)
    args.output.mkdir(parents=True, exist_ok=False)
    write_csv(args.output/'buffers.csv', buffers)
    write_csv(args.output/'groups.csv', groups)
    report = dict(capture=str(args.capture.resolve()), candidates=totals,
        scope='Sampled restriction increments; exact corners excluded from groups. '
              'Invocations are not CG iterations in multi-pass/multi-step runs.',
        model='Group size counts compressible DOFs in producer union order. '
              'Metadata model: 12 bytes per contiguous group run, 4 bytes run count, '
              '8 bytes stream lengths; exact doubles and 8-byte packet alignment. '
              'Packed payload model compares 12-bit and 16-bit fixed-width streams; '
              'it does NOT predict graph-lifting/Huffman sizes, iteration counts, '
              'total solver traffic, or MPI overhead.',
        numerical='Step ratios are ideal symmetric12 versus affine16 bounds on the '
                  'same captured increments, not actual errors or convergence guarantees. '
                  'Constant full buffers have undefined relative step ratio. '
                  'No candidate quantization is fed back to the running solver.')
    (args.output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print('group   median step ratio   p95 step ratio   packed R saving   extra metadata')
    for size, t in totals.items():
        print(size, t['symmetric12_step_ratio_p50'], t['symmetric12_step_ratio_p95'],
              t['packed_model_saving_percent'], t['metadata_extra'])
    print('These are sampled-buffer estimates, NOT measured full-codec savings.')
    print('Output:', args.output)


if __name__ == '__main__':
    main()
