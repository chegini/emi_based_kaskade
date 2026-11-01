#!/usr/bin/env python3
"""Compare five independent BDDC applications to the same first Arnoldi input."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
import re
from pathlib import Path
import shlex
import subprocess

from collect_emi_scaling_smoke import collect
from run_emi_fgmres_check import make_command as fgmres_command
from run_emi_mixed_precision_validation import digest
from run_emi_sd_cg_comparison import set_value, value

MODES = ('uncompressed', 'r16-p16', 'r8-p16', 'r16-p8', 'r8-p8')
MARKER = 'BDDC frozen probe completed: first normalized RHS; applications=1; no state advance'


def make_command(saved, mode, out):
    split = mode == 'r8-p8-split'
    if mode not in MODES and not split:
        raise ValueError('Unknown probe mode')
    cmd = fgmres_command(saved, 'uncompressed-fgmres' if mode == 'uncompressed' else 'full-fgmres16', out)
    r, p = (8, 8) if split else ((16, 16) if mode == 'uncompressed' else
                               tuple(int(x[1:]) for x in mode.split('-')))
    changes = {'--bddcPreconditionerProbe': 1, '--bddcRestrictCompressionBits': r,
               '--bddcProlongateCompressionBits': p, '--bddcCompressionBits': r,
               '--bddcRestrictZeroPreserving': int(r == 8), '--potentialDiagnostics': 0,
               '--vtk': 0, '--bddcCompressionReport': 0,
               '--bddcProlongateConstantSplit': int(split)}
    for flag, new in changes.items():
        set_value(cmd, flag, new)
    return cmd


def application_failure(text, rc):
    details = []
    for line in text.splitlines():
        line = line.strip()
        if 'what():' in line or line.startswith(('Frozen probe local error:', 'Frozen probe root error:')):
            if line not in details:
                details.append(line)
    return (f'Probe application/marker failed (exit {rc})'
            + (': '+ '; '.join(details[:5]) if details else '; inspect run.log'))


def norm(x):
    return math.sqrt(math.fsum(a*a for a in x))


def dot(x, y):
    return math.fsum(a*b for a, b in zip(x, y))


def read_probe(path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows or [int(r['dof']) for r in rows] != list(range(len(rows))):
        raise ValueError('Empty, duplicate or unordered probe DOFs')
    result = {k: [float(r[k]) for r in rows] for k in ('rhs', 'v', 'z', 'Az')}
    if not all(math.isfinite(x) for column in result.values() for x in column):
        raise ValueError('Nonfinite probe data')
    rhsnorm = norm(result['rhs'])
    if rhsnorm <= 0 or abs(norm(result['v'])-1) > 1e-12:
        raise ValueError('Invalid normalized RHS')
    if norm([v-b/rhsnorm for v, b in zip(result['v'], result['rhs'])]) > 1e-12:
        raise ValueError('Probe input is not the first normalized RHS')
    return result


def compare(data, reference):
    if data['rhs'] != reference['rhs'] or data['v'] != reference['v']:
        raise ValueError('Frozen RHS/input mismatch between independent runs')
    result = {}
    for key in ('z', 'Az'):
        denom = norm(reference[key])
        if denom == 0:
            raise ValueError('Zero reference correction/operator image')
        result[key+'_relative_error'] = norm([a-b for a, b in zip(data[key], reference[key])])/denom
    nv, naz = norm(data['v']), norm(data['Az'])
    cosine = dot(data['v'], data['Az'])/(nv*naz) if naz else None
    result.update(v_dot_z=dot(data['v'], data['z']), Az_norm=naz,
                  cosine_v_Az=cosine,
                  best_one_direction_relative_residual=(math.sqrt(max(0, 1-cosine*cosine))
                                                       if cosine is not None else 1.0))
    if any(x is not None and not math.isfinite(x) for x in result.values()):
        raise ValueError('Nonfinite comparison')
    return result


def read_traffic(text, split):
    rows = re.findall(r'^BDDC frozen probe traffic: restriction=(\d+) prolongation=(\d+) constant_header=(\d+)\s*$', text, re.M)
    if len(rows) != 1:
        raise ValueError('Missing/duplicate frozen probe packet accounting; rebuild EMI with current headers')
    r, p, header = map(int, rows[0])
    if header > p or header % 8 or bool(header) != split:
        raise ValueError('Invalid constant-header byte accounting')
    return dict(restriction_bytes=r, prolongation_bytes=p, constant_header_bytes=header,
                interface_bytes=r+p)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--reference-campaign', type=Path, required=True)
    p.add_argument('--output', type=Path)
    p.add_argument('--constant-split', action='store_true', help='add one R8/P8 constant-split candidate; no convergence run')
    p.add_argument('--reference-probe', type=Path, help='saved five-mode probe; required for --constant-split parity')
    args = p.parse_args()
    repo = args.repo.resolve(); source = args.reference_campaign.resolve()/'10cells/full-cg16'
    saved = shlex.split((source/'command.sh').read_text())
    previous = (source/'run.log').read_text(errors='replace')
    control = collect(previous, int(value(saved, '--maximumNumberOfTimeSteps')), 8, 0, 1e-7)
    if control['errors']:
        p.error('Saved CG16 control did not pass')
    model = repo/'emiModel'; sha = digest(model)
    if str(model) not in saved:
        p.error('Source executable path differs')
    modes = MODES + (('r8-p8-split',) if args.constant_split else ())
    if args.constant_split and not args.reference_probe:
        p.error('--constant-split requires --reference-probe to verify unchanged controls')
    old_vectors = {}
    if args.reference_probe:
        for mode in MODES:
            old_vectors[mode] = read_probe(args.reference_probe/mode/'bddc-frozen-probe.csv')
    for mode in modes:
        make_command(saved, mode, Path('validation-only'))
    root = (args.output or repo/'output'/('emi-frozen-probe-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+str(os.getpid()))).resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root/'environment.json').write_text(json.dumps(dict(executable_sha256=sha,
        source_command=saved, source_log_sha256=digest(source/'run.log'),
        input='first normalized RHS; fresh process/feedback/codebook state per mode',
        applications=1, restriction_passes=1, tolerance=1e-7, state_advanced=False,
        constant_split_candidate=args.constant_split,
        reference_probe=str(args.reference_probe.resolve()) if args.reference_probe else None), indent=2)+'\n')
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    env.pop('KASKADE_RESTRICTION_CAPTURE', None)
    records = {}; vectors = {}
    print('Frozen BDDC probe:', root, flush=True)
    for mode in modes:
        if digest(model) != sha:
            raise RuntimeError('Executable changed during probe')
        folder = root/mode; folder.mkdir()
        cmd = make_command(saved, mode, folder)
        (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        print('Applying '+mode+' once ...', flush=True)
        with (folder/'run.log').open('w') as f:
            rc = subprocess.run(cmd, cwd=repo, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
        text = (folder/'run.log').read_text(errors='replace')
        row = dict(status='FAIL', return_code=rc, errors=[])
        try:
            if rc != 0 or text.count(MARKER) != 1:
                raise ValueError(application_failure(text, rc))
            if 'BDDC global residual step' in text or 'BDDC FGMRES statistics step' in text or 'EMI model completed' in text:
                raise ValueError('Probe unexpectedly advanced the solver or reported a normal solve')
            data = read_probe(folder/'bddc-frozen-probe.csv')
            if len(data['v']) != control['dofs']:
                raise ValueError('Probe DOFs differ from saved control')
            if mode in old_vectors:
                if data != old_vectors[mode]:
                    raise ValueError('Saved frozen control parity failed: RHS/v/z/Az changed')
                row['saved_control_parity'] = True
            if args.constant_split and mode != 'uncompressed':
                row.update(read_traffic(text, mode == 'r8-p8-split'))
            vectors[mode] = data
            row.update(compare(data, vectors['uncompressed']))
            row['status'] = 'MEASURED'
        except (OSError, ValueError, KeyError) as e:
            row['errors'].append(str(e))
        records[mode] = row
        (root/'summary.json').write_text(json.dumps(records, indent=2, allow_nan=False)+'\n')
        print('  '+row['status'], flush=True)
        if row['status'] != 'MEASURED':
            break
    if digest(model) != sha:
        raise RuntimeError('Executable changed during probe')
    lines = ['Same first Arnoldi input; one BDDC application; no PDE step or convergence claim.',
             'mode            relative z error   relative A*z error    cos(v,A*z)    best 1-direction residual']
    for mode, row in records.items():
        if row['status'] != 'MEASURED':
            lines.append(mode+' FAIL: '+'; '.join(row['errors'])); continue
        fmt = lambda x: '--' if x is None else f'{x:.6e}'
        lines.append(f'{mode:14s} '+ '  '.join(f'{fmt(row[k]):>18s}' for k in
            ('z_relative_error', 'Az_relative_error', 'cosine_v_Az', 'best_one_direction_relative_residual')))
        if 'interface_bytes' in row:
            lines.append(f"  interface bytes={row['interface_bytes']} (R={row['restriction_bytes']}, P={row['prolongation_bytes']}); constant header included={row['constant_header_bytes']}")
    if len(records) == len(modes) and all(r['status'] == 'MEASURED' for r in records.values()):
        # Non-additivity of the two direction changes on this fixed input.
        interaction = [a-b-c+d for a,b,c,d in zip(vectors['r8-p8']['Az'], vectors['r8-p16']['Az'],
                                                vectors['r16-p8']['Az'], vectors['r16-p16']['Az'])]
        q = norm(interaction)/norm(vectors['uncompressed']['Az'])
        (root/'interaction.json').write_text(json.dumps({'Az_nonadditive_error_vs_reference_norm': q}, indent=2)+'\n')
        lines.append(f'Combined-direction non-additivity, normalized by uncompressed ||Az||: {q:.6e}')
        if args.constant_split:
            candidate = records['r8-p8-split']; baseline = records['r16-p16']
            saving = 100*(1-candidate['interface_bytes']/baseline['interface_bytes'])
            lines.append(f'First-application interface bytes saved vs16: {saving:.2f}% (not whole-solve savings)')
            lines.append('All five saved controls reproduced identical RHS/v/z/Az values.')
    lines += ['', 'Errors use uncompressed BDDC on this SAME input, not the inverse of A.',
              '8-bit restriction is exact-zero; 16-bit restriction and all prolongations use original quantizers.',
              'Separate processes retain independent fresh codebook/feedback state. Scale/minimum remain float32.',
              'This diagnoses the first application only; it does not replay later warmed-codebook states.',
              'The one-direction residual is a geometric prediction, not a solver run or convergence guarantee.',
              'Constant-split candidate, if present: quantize w*(u-c), reconstruct decoded+w*c; c is float32.',
              'New packet header: 8 bytes (float32 c + uint32 tag); setup weights reused; no extra solver message.',
              'Packet counts include metadata, padding, exact trailers and intra-rank deliveries; exclude setup/codebooks, coarse traffic and diagnostic reductions.',
              'These probes are not a convergence or timing benchmark. Previous results are untouched.']
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')
    print('\n'+'\n'.join(lines)); print('Results:', root)
    return int(len(records) != len(modes) or any(r['status'] != 'MEASURED' for r in records.values()))


if __name__ == '__main__':
    raise SystemExit(main())
