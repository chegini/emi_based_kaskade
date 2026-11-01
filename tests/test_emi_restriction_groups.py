#!/usr/bin/env python3
import tempfile
from pathlib import Path
from analyze_emi_restriction_groups import analyze, packet_cost

# Shared source DOF 1 keeps its producer group across both packets.
group = {0:0, 1:0, 2:1, 3:1}
a = packet_cost([0,1,2,3], set(), group)
b = packet_cost([0,2,1,3], set(), group)
assert a['runs'] == 2 and b['runs'] == 4
assert a['metadata_extra'] == 20 and b['metadata_extra'] == 44
assert a['packed16'] == 24 and a['packed12'] == 48
assert packet_cost([0,1,2,3,4], {4}, group)['packed12'] == 56
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    (root/'sub0-topology.csv').write_text(
        'neighbor,index,excluded\n1,0,0\n1,1,0\n2,1,0\n2,2,0\n2,3,1\n')
    (root/'sub0-call1.csv').write_text(
        'index,component,excluded,value\n0,0,0,-1\n1,0,0,1\n2,0,0,0.01\n3,0,1,1000000\n')
    buffers, groups, totals = analyze(root,[2])
    assert buffers[0]['minimum'] == -1 and buffers[0]['maximum'] == 1
    assert buffers[0]['count'] == 3
    assert groups[1]['count'] == 1 and groups[1]['span'] == 0
    assert totals[2]['packets'] == 2 and totals[2]['runs'] == 3
    # Single-value symmetric group is not an exact affine constant group.
    assert groups[1]['symmetric12_step_over_affine16'] > 0
    (root/'sub0-call1.csv').write_text(
        'index,component,excluded,value\n0,0,0,0\n1,0,0,0\n2,0,0,0\n3,0,1,1000000\n')
    _, _, totals = analyze(root,[2])
    assert totals[2]['symmetric12_step_ratio_p50'] is None
print('Restriction group analysis checks passed (ranges, exclusions, repeated groups, alignment).')
