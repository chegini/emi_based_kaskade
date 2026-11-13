#!/usr/bin/env python3
"""Check balanced sampling and compact summary arithmetic without a solver."""
import csv,importlib.util,subprocess,sys,tempfile,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('selector',HERE/'select_emi_hb_inputs.py');selector=importlib.util.module_from_spec(spec);spec.loader.exec_module(selector)
class AuditTests(unittest.TestCase):
 def test_complete_sample_and_missing_producer(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);src=root/'source';src.mkdir()
   for phase in ['restrict','prolongate']:
    for call in [1,2,3]:
     for producer in [0,1]:
      (src/f'{phase}-{call}-sub{producer}.buffer').write_text(f'KASKADE_FROZEN_V1 {phase} {call} {producer} 0\n')
   selector.select(src,root/'out',2)
   self.assertEqual(len(list((root/'out/inputs').glob('*.buffer'))),8)
   self.assertTrue(all(p.is_symlink() for p in (root/'out/inputs').iterdir()))
   (src/'restrict-2-sub1.buffer').unlink()
   with self.assertRaises(ValueError):selector.select(src,root/'bad',2)
   self.assertFalse((root/'bad').exists())
 def test_statistics_match_explicit_coefficients(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);d=root/'replay-b16-d2';d.mkdir()
   (d/'packets.csv').write_text('mode,phase,raw_bytes,wire_bytes,huffman_bytes,tail_bytes,padding_bytes,exact_bytes,plan_bytes,squared_error,squared_norm\nbefore,restrict,16,24,1,1,6,0,10,0,1\n')
   (d/'level-stats.csv').write_text('mode,phase,call,producer,peer,level,count,minimum,maximum,sum_squares,tail_bits\nbefore,restrict,1,0,1,2,2,-2,4,20,3\nbefore,restrict,2,0,1,2,1,3,3,9,2\n')
   subprocess.run([sys.executable,str(HERE/'summarize_emi_hb_audit.py'),str(root)],check=True,stdout=subprocess.DEVNULL)
   with (root/'coefficient-levels.csv').open() as f:r=next(csv.DictReader(f))
   self.assertEqual(int(r['count']),3);self.assertEqual(float(r['minimum']),-2);self.assertEqual(float(r['maximum']),4)
   self.assertAlmostEqual(float(r['rms']),(29/3)**.5);self.assertEqual(int(r['tail_bits']),5)
if __name__=='__main__':unittest.main()
