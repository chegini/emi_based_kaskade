#!/usr/bin/env python3
"""Configuration, route weighting and environment separation regressions."""
import os,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
import run_emi_message_entropy as e
from check_emi_entropy_observer import fixture_source
class Workflow(unittest.TestCase):
 def test_actual_production_fixture(self):
  source=Path(__file__).with_name('test_bddc_split_recenter.cpp').read_text()
  result=fixture_source(source)
  selector='\n   if(std::getenv("EMI_ENTROPY_TEST_A")){t.setRestrictTransformAfterQuantization(true);t.setProlongateTransformAfterQuantization(true);}'
  self.assertEqual(result.count(selector),1)
  self.assertEqual(result.replace(selector,''),source)
 def test_fixture_formatting(self):
  for text in ['   t.setProlongateZigzagEncoding(true);',
               '   t.setRestrictZigzagEncoding(true);t.setProlongateZigzagEncoding(true);',
               't . setProlongateZigzagEncoding ( true ) ;']:
   with self.subTest(text=text):
    self.assertIn('EMI_ENTROPY_TEST_A',fixture_source(text))
  for text in ['', 't.setProlongateZigzagEncoding(true);t.setProlongateZigzagEncoding(true);']:
   with self.assertRaises(ValueError):fixture_source(text)
 def test_configs(self):
  a=SimpleNamespace(bits=[8,12,16,24,32],depths=[1,2,3,4],depth_split_study=True,
    case='pepe',ranks=16,nodes=8,steps=3,case_minutes=60,
    slurm_container=Path('/fake/image.sif'),kaskade=Path('/fake/overlay'),
    executable=Path('/fake/root/build/emiModel'),output=Path('/fake/root/results'))
  cases=e.h.configurations(a);self.assertEqual(len(cases),121)
  for case in cases:
   cfg,cmd=e.command(a,case,a.output/case['name'],Path('/fake/original'))
   self.assertIn('/fake/original:/fake/original',cmd)
   self.assertEqual(cmd[cmd.index('--masterWorker')+1],'0')
   self.assertEqual(cmd[cmd.index('--bddcTransformFirstRecentering')+1],str(int(case['variant']=='B-recenter')))
   self.assertEqual(cmd[cmd.index('--bddcGraphLiftingAfterQuantization')+1],str(int(case['variant']=='A')))
   self.assertEqual(cmd[cmd.index('--bddcProlongateConstantSplit')+1],str(int(case['split']=='on')))
   self.assertEqual(cmd[cmd.index('--bddcHistogram')+1],'0')
   if case['bits']:
    self.assertEqual(cmd[cmd.index('--bddcCompressionBits')+1],str(case['bits']))
    self.assertEqual(cmd[cmd.index('--bddcGraphLiftingLevels')+1],str(case['depth']))
   self.assertNotIn('--bddcHistogramSteps',cmd)
   env=e.environment(case,a.output/case['name']);self.assertEqual(env['KASKADE_TRANSFER_TIMING'],'0')
   self.assertEqual(env['EMI_MESSAGE_ENTROPY_PREFIX'],env['APPTAINERENV_EMI_MESSAGE_ENTROPY_PREFIX'])
 def test_route_weighting(self):
  rows=[]
  for phase in ('restrict','prolongate'):
   for stage in ('quantized_input','encoder_coefficients','huffman_bitlength','routed_bytes'):
    for route,n,h in [('same_rank',100,9),('same_node',3,1),('inter_node',1,3)]:
     rows.append(dict(phase=phase,stage=stage,route=route,symbols=n,messages=1,ideal_message_entropy_bits=n*h))
  out=e.aggregate(dict(name='B',variant='B',bits=16,depth=1),rows)
  self.assertEqual(len(out),4)
  self.assertAlmostEqual(out[0]['encoder_H0_bits_per_coefficient'],1.5)
  self.assertAlmostEqual(out[1]['encoder_H0_bits_per_coefficient'],3)
  self.assertEqual(out[0]['quantized_coefficients'],4)
  self.assertEqual(out[0]['effective_routed_bits_per_quantized_coefficient'],8)
if __name__=='__main__':unittest.main()
