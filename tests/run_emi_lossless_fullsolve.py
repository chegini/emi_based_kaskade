#!/usr/bin/env python3
"""Matched one-time-step FCG solves: unchanged numerics, measured packet bytes."""
import argparse
import copy
import csv
import json
import math
import os
from pathlib import Path
import re
import shlex

import run_emi_solver_study as study

FILES = ('potential-state-step0.csv', 'potential-state-step1.csv',
         'bddc-convergence-step1.csv')


def compare(old, new):
    for filename in FILES:
        if (old['folder']/filename).read_bytes() != (new['folder']/filename).read_bytes():
            raise ValueError(filename+' changed with packing')
    for key in ('iterations_per_step', 'relative_residuals', 'original_bytes',
                'residual_lines', 'routes', 'raw_sent_bytes', 'transmitted_sent_bytes'):
        if old[key] != new[key]:
            raise ValueError(key+' changed with packing')
    if old['actual_bytes'] != old['original_bytes']:
        raise ValueError('Baseline transport counters differ')
    if old['original_bytes'] <= 0:
        raise ValueError('No inter-rank packets exercised')
    return 100*(1-new['actual_bytes']/old['actual_bytes'])


def validate_trace(samples):
    expected = {'time_step','iteration','step_length','residual_norm','rhs_norm','relative_residual','target'}
    if not samples:
        raise ValueError('Empty convergence trace')
    for row in samples:
        if set(row) != expected:
            raise ValueError('Unexpected convergence trace columns')
        for key, value in row.items():
            # Production CSV leaves alpha absent at the initial residual.
            if key == 'step_length' and row['iteration'] == '0' and value == '':
                continue
            if not math.isfinite(float(value)):
                raise ValueError('Nonfinite convergence trace')
    if [int(row['iteration']) for row in samples] != list(range(len(samples))):
        raise ValueError('Missing or unordered convergence samples')


def inspect(folder, cfg, rc, timed_out):
    log = (folder/'run.log').read_text(errors='replace')
    result = study.collect(log, 1, cfg.ranks, rc, 1e-7)
    if timed_out or result['errors']:
        raise ValueError(str(folder)+': '+str(result['errors'])+' timeout='+str(timed_out))
    result['routes'] = study.routing_counts(log, cfg.expected_nodes)
    if sum(r['bytes'] for r in result['routes']) != result['transmitted_sent_bytes']:
        raise ValueError('Logical routing counters disagree')
    matches = re.findall(r'^BDDC packet transport bytes: original=(\d+) actual=(\d+) packing=([01])$', log, re.M)
    if len(matches) != 1 or int(matches[0][2]) != cfg.lossless_packing:
        raise ValueError('Missing, duplicate or incorrect packet transport summary')
    result['original_bytes'], result['actual_bytes'] = map(int, matches[0][:2])
    if cfg.lossless_packing:
        length_summary = re.findall(r'^BDDC packet length compaction bytes: wide=(\d+) actual=(\d+)$', log, re.M)
        if len(length_summary) != 1:
            raise ValueError('Missing or duplicate length compaction summary')
        wide, actual = map(int, length_summary[0])
        if actual != result['actual_bytes'] or actual > wide:
            raise ValueError('Invalid length compaction counters')
        result['wide_bytes'] = wide
        tail_summary = re.findall(r'^BDDC packet tail omission bytes: two_lengths=(\d+) actual=(\d+)$', log, re.M)
        if len(tail_summary) != 1:
            raise ValueError('Missing or duplicate tail omission summary')
        two, tail_actual = map(int, tail_summary[0])
        if tail_actual != actual or not actual <= two <= wide:
            raise ValueError('Invalid tail omission counters')
        result['two_length_bytes'] = two
        plans = re.findall(r'^BDDC packet endpoint plan bytes: explicit=(\d+) actual=(\d+) plan_ids=(\d+) omitted_ids=(\d+) plan_messages=(\d+) reused_messages=(\d+)$', log, re.M)
        if len(plans) != 1:
            raise ValueError('Missing/duplicate endpoint plan summary')
        explicit, plan_actual, ids, omitted, messages, reused = map(int, plans[0])
        count_summary = re.findall(r'^BDDC packet exact count omission bytes: with_counts=(\d+) actual=(\d+) saved=(\d+)$', log, re.M)
        if len(count_summary) != 1:
            raise ValueError('Missing/duplicate exact-count summary')
        with_counts, count_actual, count_saved = map(int, count_summary[0])
        small_summary = re.findall(r'^BDDC packet small metadata bytes: previous=(\d+) actual=(\d+) saved=(\d+) messages=(\d+)$', log, re.M)
        if len(small_summary) != 1:
            raise ValueError('Missing/duplicate small-metadata summary')
        previous_small, small_actual, small_saved, small_messages = map(int, small_summary[0])
        if small_actual != actual or previous_small-actual != small_saved or small_saved < 0:
            raise ValueError('Invalid small-metadata accounting')
        dict_count_summary = re.findall(r'^BDDC packet dictionary count bytes: wide=(\d+) actual=(\d+) saved=(\d+) messages=(\d+)$', log, re.M)
        if len(dict_count_summary) != 1: raise ValueError('Missing/duplicate dictionary count summary')
        count_wide, dc_actual, dc_saved, dc_messages = map(int, dict_count_summary[0])
        if dc_actual != actual or count_wide < actual or count_wide-previous_small != dc_saved or dc_saved != 3*dc_messages:
            raise ValueError('Invalid dictionary count accounting')
        dictionary_summary = re.findall(r'^BDDC packet dictionary reference bytes: wide=(\d+) actual=(\d+) saved=(\d+) messages=(\d+) byte_refs=(\d+) wide_refs=(\d+)$', log, re.M)
        if len(dictionary_summary) != 1: raise ValueError('Missing/duplicate dictionary reference summary')
        dict_wide, dict_actual, dict_saved, dict_messages, byte_refs, wide_refs = map(int, dictionary_summary[0])
        if dict_actual != actual or not actual <= count_wide <= dict_wide or dict_wide-count_wide != dict_saved or dict_saved != 3*byte_refs or bool(dict_messages) != bool(byte_refs+wide_refs) or dc_messages > dict_messages or dc_messages > byte_refs or bool(dc_messages) != bool(byte_refs):
            raise ValueError('Invalid dictionary reference accounting')
        record_summary = re.findall(r'^BDDC packet record length bytes: wide=(\d+) actual=(\d+) saved=(\d+) records=(\d+) shortened=(\d+) le255=(\d+) le4095=(\d+) le32767=(\d+) larger=(\d+)$', log, re.M)
        if len(record_summary) != 1:
            raise ValueError('Missing/duplicate record length summary')
        record_wide, record_actual, record_saved, records, shortened, tiny, small, medium, large = map(int, record_summary[0])
        if (record_actual != actual or not actual <= dict_wide < record_wide <= with_counts or
                record_wide-dict_wide != record_saved or record_saved != 2*shortened or
                records != tiny+small+medium+large or not 0 < shortened <= tiny+small+medium):
            raise ValueError('Invalid or unexercised short record lengths')
        if (count_actual != actual or not actual < with_counts <= explicit or
                with_counts-record_wide != count_saved or count_saved % 4):
            raise ValueError('Invalid or unexercised exact-count omission')
        if (plan_actual != actual or not actual <= explicit <= two or
                explicit-with_counts != omitted or ids % 8 or omitted % 8 or
                not all((ids, omitted, messages, reused))):
            raise ValueError('Invalid or unexercised endpoint plans')
        result['endpoint'] = dict(explicit_id_packer_bytes=explicit,
            endpoint_plan_id_bytes=ids, endpoint_ids_bytes_saved=omitted,
            endpoint_ids_savings_percent=100*omitted/explicit,
            endpoint_plan_messages=messages, endpoint_plan_reused_messages=reused,
            with_exact_count_packer_bytes=with_counts, exact_count_bytes_saved=count_saved,
            exact_count_savings_percent=100*count_saved/with_counts,
            with_wide_record_length_bytes=record_wide,record_length_bytes_saved=record_saved,
            record_length_savings_percent=100*record_saved/record_wide,
            compact_records=records,shortened_records=shortened,
            record_bytes_0_255=tiny,record_bytes_256_4095=small,
            record_bytes_4096_32767=medium,record_bytes_over_32767=large,
            with_wide_dictionary_reference_bytes=dict_wide,dictionary_reference_bytes_saved=dict_saved,
            dictionary_reference_savings_percent=100*dict_saved/dict_wide if dict_wide else 0,
            dictionary_messages=dict_messages,byte_dictionary_references=byte_refs,
            wide_dictionary_references=wide_refs,
            with_wide_dictionary_count_bytes=count_wide,dictionary_count_bytes_saved=dc_saved,
            dictionary_count_savings_percent=100*dc_saved/count_wide if count_wide else 0,
            short_dictionary_counts=dc_messages,
            before_small_metadata_bytes=previous_small,
            small_metadata_bytes_saved=small_saved,
            small_metadata_saving_percent=100*small_saved/previous_small if previous_small else 0,
            small_metadata_messages=small_messages)

    for step in (0, 1):
        study.read(folder/f'potential-state-step{step}.csv')
    with (folder/FILES[2]).open() as stream:
        samples = list(csv.DictReader(stream))
    validate_trace(samples)
    if int(samples[-1]['iteration']) != result['iterations_per_step'][0]:
        raise ValueError('Trace iteration count differs from summary')
    result['residual_lines'] = [s for s in log.splitlines() if s.startswith('BDDC global residual step ')]
    result['folder'] = folder
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--executable', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--kaskade', type=Path, required=True)
    p.add_argument('--slurm-container', type=Path, required=True)
    p.add_argument('--case', default='weak_32_h16')
    p.add_argument('--ranks', type=int, default=4)
    p.add_argument('--expected-nodes', type=int, default=2)
    p.add_argument('--bits', type=int, nargs='+', choices=[8,12,16,24,32], default=[16])
    p.add_argument('--lifting-levels', type=int, choices=range(65), default=1)
    p.add_argument('--split', choices=['off','on'], default='off')
    p.add_argument('--case-minutes', type=int, default=30)
    a = p.parse_args()
    cfg = study.arguments(['--repo', str(a.repo), '--executable', str(a.executable),
        '--kaskade', str(a.kaskade), '--slurm-container', str(a.slurm_container),
        '--case', a.case, '--ranks', str(a.ranks), '--expected-nodes', str(a.expected_nodes),
        '--bits', *map(str,a.bits), '--solvers', 'fcg', '--split', a.split,
        '--case-minutes', str(a.case_minutes), '--lifting-levels', str(a.lifting_levels)])
    if not os.environ.get('SLURM_JOB_ID'):
        p.error('Run inside a Slurm allocation')
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    sources = [cfg.executable, cfg.repo/'EMI_model.cpp', *cfg.repo.glob('*.hh'),
               Path(__file__), cfg.repo/'tests/run_emi_solver_study.py',
               cfg.repo/'tests/run_emi_histogram_start_check.py',
               cfg.repo/'tests/run_emi_mpi_routing_check.py']
    sources += [cfg.kaskade/'mg'/s for s in ('bddc.hh','bddc.hpp','bddcSpaceTransfer.hh',
                 'bddcSpaceTransfer.hpp','bddcPacketPacking.hh')]
    sources += [Path(s) for s in study.mesh_inputs(cfg.repo, cfg.case).values()]
    hashes = {str(s):study.digest(s) for s in sources}
    (root/'manifest.json').write_text(json.dumps(hashes,indent=2)+'\n')
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for key in ('KASKADE_FROZEN_TRANSFER_DIR','KASKADE_TRANSFER_AUDIT_DIR'):
        env.pop(key, None)
    rows = []
    for bits in a.bits:
        for order in ('before','after'):
            pair = []
            for enabled in (0,1):
                run_cfg = copy.copy(cfg)
                run_cfg.lossless_packing = enabled
                run_cfg.transform_order = order
                name = f'fcg-{order}-{bits}-packing-{enabled}'
                folder = root/name
                folder.mkdir()
                case = dict(name=name, solver='fcg', bits=bits, split=a.split, direct=False)
                cmd = study.make_command(run_cfg, case, folder, ['mpirun'])
                (folder/'command.sh').write_text(shlex.join(cmd)+'\n')
                print('Running '+name, flush=True)
                rc, timed_out, _ = study.execute(cmd, folder, run_cfg, env)
                pair.append(inspect(folder, run_cfg, rc, timed_out))
            old,new = pair
            saving = compare(old,new)
            row = dict(bits=bits, order=order, iterations_off=old['iterations_per_step'][0],
                iterations_on=new['iterations_per_step'][0], true_residual=new['relative_residuals'][0],
                original_routed_bytes=old['actual_bytes'], packed_routed_bytes=new['actual_bytes'],
                additional_savings_percent=saving, parity='EXACT_SOLUTION_TRACE_AND_RESIDUALS')
            wide = new['wide_bytes']
            two = new['two_length_bytes']
            explicit = new['endpoint']['explicit_id_packer_bytes']
            row.update(previous_packer_bytes=wide,
                stream_length_bytes_saved=wide-two,
                stream_length_savings_percent=100*(1-two/wide) if wide else 0,
                two_length_packer_bytes=two,tail_length_bytes_saved=two-explicit,
                tail_length_savings_percent=100*(1-explicit/two) if two else 0)
            row.update(new['endpoint'])
            rows.append(row)
            with (root/'summary.csv').open('w') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(row))
                writer.writeheader()
                writer.writerows(rows)
            print('PASS: '+json.dumps(row), flush=True)
    if hashes != {str(s):study.digest(s) for s in sources}:
        raise ValueError('Sources, input data or executable changed during validation')
    (root/'PASS').write_text('All matched EMI FCG full-solve comparisons passed.\n')
    print('Validation complete: '+str(root), flush=True)


if __name__ == '__main__':
    main()
