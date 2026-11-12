import csv
import json
from pathlib import Path
import tempfile
import unittest
from compare_emi_multilevel_hb import summarize

class Report(unittest.TestCase):
    def test_worse_compression_and_failed_cases_remain_visible(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for depth in (1,0):
                for order in ('before','after'):
                    folder=root/f'levels-{depth}'/order
                    folder.mkdir(parents=True)
                    runs={'fcg-uncompressed':dict(status='PASS',packet_transport=dict(actual_bytes=1000))}
                    for bits in (8,12,16,24,32):
                        name=f'fcg-{order}-{bits}-split-on'
                        sub=folder/name;sub.mkdir()
                        sub.joinpath('run.log').write_text('BDDC integer lifting setup rank=0 local_plan_bytes=10 exported_remote_bytes=20\n')
                        runs[name]=dict(status='FAIL' if depth==0 and bits==8 else 'PASS',
                            packet_transport=dict(actual_bytes=200 if depth==1 else 210),
                            iterations_per_step=[10],error_vs_direct=dict(relative_l2=1e-6))
                    folder.joinpath('summary.json').write_text(json.dumps(dict(runs=runs)))
            summarize(root)
            with (root/'hb-comparison.csv').open() as stream:
                rows=list(csv.DictReader(stream))
            bad=next(r for r in rows if r['levels']=='0' and r['bits']=='8')
            self.assertEqual(bad['actual_routed_bytes'],'')
            self.assertEqual(bad['bytes_saved_vs_single'],'')
            worse=next(r for r in rows if r['levels']=='0' and r['bits']=='12')
            self.assertEqual(worse['bytes_saved_vs_single'],'-10')
            self.assertAlmostEqual(float(worse['saving_gain_percentage_points']),-1)
            self.assertEqual(worse['exported_setup_bytes'],'20')

if __name__=='__main__': unittest.main()
