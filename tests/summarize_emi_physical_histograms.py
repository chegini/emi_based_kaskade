#!/usr/bin/env python3
"""Export first-call EMI physical histograms; no spatial-frequency inference."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path


def summarize(root):
    data=root/'physical-histograms'
    data.mkdir(exist_ok=True)
    result={}
    solver_label=None
    for onset in (1,6):
        campaign=root/f'onset-{onset}'
        environment=json.loads((campaign/'environment.json').read_text())
        solver=environment['solver']
        label={'fcg':'FCG', 'cg':'CG', 'fgmres':'FGMRES(30)'}[solver]
        if solver_label is not None and solver_label!=label:
            raise ValueError('Mixed solver families across compression onsets')
        solver_label=label
        check=json.loads((campaign/'summary.json').read_text())
        if check['errors'] or check['checks'].get('observer')!='EXACT_SNAPSHOTS_TRACES_ITERATIONS_AND_BYTES':
            raise ValueError('Observer parity gate failed: '+str(campaign))
        base=campaign/'scheduled-histograms'/'symbol-histograms'/f'step{onset}'
        for phase in ('restrict','prolongate'):
            values={'input':[], 'coefficient':[]}
            files=sorted(base.glob(f'subdomain*/physical/{phase}-exchange1-neighbor*.csv'))
            if not files: raise ValueError('Missing physical capture '+str(base)+' '+phase)
            provenance={}
            for path in files:
                raw=path.read_bytes();provenance[str(path.relative_to(root))]=hashlib.sha256(raw).hexdigest()
                rows=list(csv.DictReader(raw.decode().splitlines()))
                symbols=Counter(int(r['quantized']) for r in rows)
                original=path.parent.parent/path.name
                with original.open() as f:
                    expected=Counter({int(r['symbol']):int(r['count']) for r in csv.DictReader(f)
                                      if r['stage']=='quantized_input'})
                if symbols!=expected: raise ValueError('Physical/code populations disagree: '+str(path))
                for key in values: values[key].extend(float(r[key]) for r in rows)
            if not values['input'] or not all(math.isfinite(v) for xs in values.values() for v in xs):
                raise ValueError('Empty or nonfinite population')
            maximum=max(abs(v) for xs in values.values() for v in xs)
            # All-zero data use a plotting scale of one; actual maximum remains recorded as zero.
            scale=maximum or 1.0
            stats={'solver':solver,'n':len(values['input']),'directed_packets':len(files),'maximum':maximum,
                   'plot_scale':scale,'source_sha256':provenance}
            for stage,xs in values.items():
                stats[stage]={'min':min(xs),'max':max(xs),
                             'fraction_abs_le_1pct_scale':sum(abs(v)<=.01*scale for v in xs)/len(xs)}
                for kind in ('signed','magnitude'):
                    lo=-1 if kind=='signed' else 0
                    width=(1-lo)/64
                    counts=[0]*64
                    for v in xs:
                        x=v/scale if kind=='signed' else abs(v)/scale
                        counts[min(63,max(0,int((x-lo)/width)))]+=1
                    assert sum(counts)==len(xs)
                    with (data/f'step{onset}-{phase}-{stage}-{kind}.csv').open('w') as f:
                        w=csv.writer(f);w.writerow(['x','count','percent'])
                        w.writerows((lo+(i+.5)*width,n,100*n/len(xs)) for i,n in enumerate(counts))
            result[f'step{onset}-{phase}']=stats
    (data/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (data/'README.txt').write_text(f'Pepe {solver_label}, 16 bits, transform first, splitting ON.\n'
        'First transfer call per producer at each separately selected compression onset.\n'
        'Input is weighted and, on prolongation, already centered by constant splitting.\n'
        'Coefficient is the actual floating quantizer input. Exact corners are excluded.\n'
        'Repeated directed packet copies are counted; coarse/fine modes are not classified.\n'
        'Each phase/onset uses a common scale for input and coefficients; 64 equal-width bins.\n'
        'Histogram occurrence frequency is not spatial frequency. This is not a runtime benchmark.\n')
    template=Path(__file__).with_name('emi_physical_histograms.tex')
    (data/'histograms.tex').write_text(template.read_text().replace('@@SOLVER@@',solver_label))
    print('Verified and exported physical histogram data:',data)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    summarize(p.parse_args().output.resolve())
