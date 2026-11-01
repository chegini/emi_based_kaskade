#!/usr/bin/env python3
"""Test the scientific acceptance criteria, without launching an EMI simulation."""
import copy
import unittest
from run_emi_zero_preserving_check import assess

class Acceptance(unittest.TestCase):
    def test_acceptance(self):
        baseline={'status':'PASS','restriction_passes':1,'relative_residual':7e-8,
                  'logical_transferred_bytes':2430000,
                  'comparison_to_direct':{'relative L2 error':2.8e-6,'max absolute error':3e-7}}
        candidate=copy.deepcopy(baseline)
        candidate['logical_transferred_bytes']=2000000
        self.assertEqual(assess(baseline,candidate),[])
        for key,value in [('restriction_passes',2),('relative_residual',1.01e-7),
                          ('status','FAIL'),('logical_transferred_bytes',2430000),
                          ('logical_transferred_bytes',None)]:
            bad=copy.deepcopy(candidate);bad[key]=value
            self.assertTrue(assess(baseline,bad),key)
        for key in ['relative L2 error','max absolute error']:
            bad=copy.deepcopy(candidate);bad['comparison_to_direct'][key]*=2
            self.assertTrue(assess(baseline,bad),key)
        bad=copy.deepcopy(candidate);bad['relative_residual']=float('nan')
        self.assertTrue(assess(baseline,bad))
        # Both the historical and the matched baseline must be beaten.
        smaller=copy.deepcopy(baseline);smaller['logical_transferred_bytes']=1900000
        self.assertTrue(assess(smaller,candidate))

if __name__=='__main__':unittest.main()
