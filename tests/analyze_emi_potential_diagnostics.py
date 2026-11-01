#!/usr/bin/env python3
"""Read full-precision DOF snapshots; never adjust solver states or acceptance criteria."""
import argparse
import csv
import json
import math
from pathlib import Path


def read(path):
    with path.open() as stream:
        rows=list(csv.DictReader(stream))
    if not rows or [int(r['dof']) for r in rows]!=list(range(len(rows))):
        raise ValueError('Missing or unordered global DOFs: '+str(path))
    data={k:[float(r[k]) for r in rows] for k in rows[0] if k!='dof'}
    if any(not math.isfinite(x) for col in data.values() for x in col):
        raise ValueError('Nonfinite snapshot: '+str(path))
    return data


def dot(x,y):
    if len(x)!=len(y): raise ValueError('DOF counts differ')
    return math.fsum(a*b for a,b in zip(x,y))


def error_metrics(u,reference):
    if len(u)!=len(reference): raise ValueError('DOF counts differ')
    e=[a-b for a,b in zip(u,reference)]; n=len(e)
    mean=math.fsum(e)/n
    square=dot(e,e); centered=math.fsum((x-mean)**2 for x in e)
    ref=dot(reference,reference)
    return dict(mean_error=mean,absolute_l2=math.sqrt(square),
                relative_l2=math.sqrt(square/ref) if ref else None,
                centered_relative_l2=math.sqrt(centered/ref) if ref else None,
                max_absolute_error=max(map(abs,e)),
                mean_fraction_of_squared_error=n*mean*mean/square if square else 0.0)


def write_csv(path,rows):
    if not rows: return
    with path.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def analyze(root):
    env=json.loads((root/'environment.json').read_text())
    summary=json.loads((root/'summary.json').read_text())
    out=root/'potential-analysis'
    out.mkdir(exist_ok=False)
    all_rows=[]; operators=[]; missing=[]
    for case,record in summary.items():
        base=root/case
        direct_op=read(base/'direct/potential-operator.csv')
        weights=direct_op['a_transpose_one']; n=len(weights)
        anchor=math.fsum(weights)
        for mode in record['runs']:
            folder=base/mode
            opfile=folder/'potential-operator.csv'
            if not opfile.exists():
                missing.append(str(opfile));continue
            op=read(opfile)
            if op!=direct_op:
                raise ValueError('Operator diagnostic differs from direct: '+str(opfile))
            expected=record['runs'][mode].get('metrics',{}).get('dofs')
            if expected is not None and n!=expected:
                raise ValueError('Snapshot DOF count differs from solver log: '+str(opfile))
            operators.append(dict(case=case,mode=mode,dofs=n,
                one_A_one=anchor,constant_rayleigh=anchor/n,
                norm_A_one=math.sqrt(dot(op['a_one'],op['a_one'])),
                mean_diagonal=math.fsum(op['diagonal'])/n,
                note='Constant-direction diagnostic, not an eigenvalue or condition estimate'))
            previous=None; previous_ref=None
            for step in range(env['steps']+1):
                path=folder/f'potential-state-step{step}.csv'
                refpath=base/'direct'/path.name
                if not path.exists() or not refpath.exists():
                    missing.append(str(path if not path.exists() else refpath));break
                state=read(path);ref=read(refpath)
                u=state['u']; reference=ref['u']
                if step==0 and u!=reference:
                    raise ValueError('Initial states do not match direct: '+str(path))
                errors=error_metrics(u,reference)
                one_A_du=dot(weights,state['correction'])
                sum_rhs=math.fsum(state['rhs']); sum_r=math.fsum(state['residual'])
                one_A_u=dot(weights,u)
                row=dict(case=case,mode=mode,step=step,time=step*env['physics']['dt'],
                         **errors,
                         mean_u=math.fsum(u)/n,mean_direct=math.fsum(reference)/n,
                         mean_correction_error=math.fsum(a-b for a,b in zip(state['correction'],ref['correction']))/n,
                         mean_error_increment=(errors['mean_error']-
                             math.fsum(a-b for a,b in zip(previous,previous_ref))/n) if previous is not None else None,
                         one_A_u=one_A_u,
                         anchor_weighted_error_proxy=dot(weights,[a-b for a,b in zip(u,reference)])/anchor if anchor else None,
                         sum_rhs=sum_rhs,sum_residual=sum_r,
                         sum_residual_direct=math.fsum(ref['residual']),
                         constant_balance_closure=sum_rhs-one_A_du-sum_r,
                         rhs_plus_previous_anchor=sum_rhs+dot(weights,previous) if previous is not None else None,
                         post_anchor_plus_residual=one_A_u+sum_r if step else None)
                all_rows.append(row)
                previous=u;previous_ref=reference
    write_csv(out/'per-step-errors.csv',all_rows)
    write_csv(out/'constant-mode.csv',operators)
    report=dict(complete=not missing,missing=missing,
        interpretation=[
            'Arithmetic means and Euclidean norms of unique global algebraic DOFs, not volume-weighted FE norms or nonconforming VTK samples.',
            'Centered errors are diagnostic only; raw acceptance criteria and all states are unchanged.',
            'Operator column sums define the anchor-weighted proxy; interpreting it as a boundary average requires the constant-mode cancellations.',
            'Step zero is the initial state with zero diagnostic RHS/correction, not a solve.',
            'Snapshots are after accepted updates; the RHS of each run depends on its own preceding state.',
            'No membrane-voltage or gradient error is inferred from a small centered potential error.',
            'Observe constant_balance_closure to assess summation roundoff; rhs_plus_previous_anchor tests the expected boundary balance.'])
    (out/'interpretation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Full-precision potential diagnostics:',out)
    print('case mode step mean-error raw-relative-L2 centered-relative-L2')
    for row in all_rows:
        if row['mode']=='direct': continue
        print(row['case'],row['mode'],row['step'],
              *(f'{row[k]:.6e}' if row[k] is not None else '--'
                for k in ['mean_error','relative_l2','centered_relative_l2']))
    if missing: print('Incomplete snapshots; see interpretation.json')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    args=p.parse_args()
    analyze(args.run.resolve())
