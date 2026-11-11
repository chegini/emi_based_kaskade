#!/usr/bin/env python3
"""Compare new split-FCG runs with the rounded values printed in Table 5.

Savings are inter-node logical R/P payload, NOT packed routed-message bytes.
Previous errors are printed/rounded values, so differences are approximate.
"""
import csv
from pathlib import Path
import re
import sys


def compare(root):
    previous = list(csv.DictReader((root/'previous.csv').open()))
    result = []
    for order, index in [('before', 0), ('after', 1)]:
        current = {r['case']: r for r in csv.DictReader((root/order/'results.csv').open())}
        for old in previous:
            bits = old['bits']
            new = current[f'fcg-{order}-{bits}-split-on']
            match = re.search(r'([\d.]+)\\times10\^\{(-?\d+)\}', old['e'+order])
            old_error = float(match[1])*10**int(match[2]) if match else None
            old_iter = int(old['kpair'].split('/')[index].strip())
            old_saved = float(old['s'+order])
            valid = new['status'] == 'PASS'
            result.append(dict(bits=bits, order=order, status=new['status'],
                old_iterations=old_iter, new_iterations=new['iterations'],
                iteration_change=int(new['iterations'])-old_iter if valid else '',
                old_relative_l2_rounded=old_error, new_relative_l2=new['raw_relative_l2_direct'],
                old_inter_node_saved_percent=old_saved,
                new_inter_node_saved_percent=new['inter_node_bytes_saved_percent'],
                saving_change_percentage_points=float(new['inter_node_bytes_saved_percent'])-old_saved if valid else '',
                new_inter_node_bytes=new['inter_node_bytes'],
                new_wrapper_seconds=new['solver_wrapper_seconds']))
    with (root/'comparison.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(result[0]))
        writer.writeheader(); writer.writerows(result)
    print('Table 5 comparison (old values rounded; inter-node logical payload):')
    for r in result:
        print(f"{r['order']:6} {r['bits']:>2} bits {r['status']:4}: iterations {r['old_iterations']} -> {r['new_iterations']}; savings {r['old_inter_node_saved_percent']} -> {r['new_inter_node_saved_percent']} %")
    print('Comparison:', root/'comparison.csv')


if __name__ == '__main__':
    compare(Path(sys.argv[1]))
