"""Independent byte-count/entropy validation of the production observer export."""
from collections import Counter
import csv
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

from test_emi_compression_report import kaskade_include
import run_emi_partition_pilot as pilot


class SerializedBytes(unittest.TestCase):
    def test_production_observer(self):
        repo=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory(prefix='emi-byte-entropy-') as tmp:
            root=Path(tmp);exe=root/'test';out=root/'report'
            subprocess.run(shlex.split(os.environ.get('CXX','c++'))+[
                '-std=c++20','-O2','-Wall','-Wextra','-Werror','-I',kaskade_include(repo),
                str(repo/'tests/test_bddc_serialized_byte_entropy.cpp'),'-o',str(exe)],check=True)
            result=subprocess.run([str(exe),str(out)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            with (out/'serialized-message-byte-entropy.csv').open() as f:rows=list(csv.DictReader(f))
            with (out/'serialized-message-byte-counts.csv').open() as f:bins=list(csv.DictReader(f))
            self.assertEqual(len(rows),6);self.assertEqual(len(bins),6*256)
            expected={('restriction','same_node'):Counter({0:5}),
                      ('restriction','inter_node'):Counter({255:5}),
                      ('prolongation','same_node'):Counter([0,0,0,1,0,0,0,2,16,0,0,0,0,0,0,0,
                                                           0,255,127,128,42,42,42,0]),
                      ('prolongation','inter_node'):Counter()}
            for phase in ['restriction','prolongation']:
                expected[phase,'inter_rank']=expected[phase,'same_node']+expected[phase,'inter_node']
            for r in rows:
                phase,scope=r['phase'],r['scope'];counts=expected[phase,scope];n=sum(counts.values())
                self.assertEqual(int(r['total_message_bytes']),n)
                actual={int(b['byte_value']):int(b['count']) for b in bins if b['phase']==phase and b['scope']==scope}
                self.assertEqual(actual,{b:counts[b] for b in range(256)})
                self.assertEqual(r['entropy_aggregation'],'pooled_all_messages_all_ranks')
                if n:
                    h=-math.fsum(v/n*math.log2(v/n) for v in counts.values())
                    self.assertAlmostEqual(float(r['serialized_byte_entropy_bits_per_byte']),h,places=13)
                else:
                    self.assertEqual(r['serialized_byte_entropy_bits_per_byte'],'')
                    self.assertEqual(r['mean_message_bytes'],'');self.assertEqual(r['metric_status'],'EMPTY')
            self.assertFalse((out/'disabled').exists())
            transport=[]
            for (phase,scope),counts in expected.items():
                messages=(2 if scope=='inter_rank' else 1) if phase=='restriction' else (0 if scope=='inter_node' else 1)
                transport.append(dict(phase=phase,scope=scope,sent_bytes=sum(counts.values()),mpi_messages=messages))
            values=pilot.serialized_entropy(out,transport)
            self.assertEqual(values['restriction_inter_rank_byte_entropy_bits_per_byte'],1)
            self.assertIsNone(values['prolongation_inter_node_byte_entropy_bits_per_byte'])
            transport[0]['sent_bytes']+=1
            with self.assertRaisesRegex(ValueError,'population differs'):
                pilot.serialized_entropy(out,transport)
            print('PASS: synthetic observer math/count/export and actual-population checks (not EMI measurements)')


if __name__=='__main__':unittest.main()
