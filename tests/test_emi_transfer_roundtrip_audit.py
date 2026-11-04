import importlib.util
from pathlib import Path
import tempfile
import unittest
spec=importlib.util.spec_from_file_location('audit',Path(__file__).with_name('run_emi_transfer_roundtrip_audit.py'))
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

class PacketAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        for phase in ('restrict','prolongate'):
            for stage in audit.STAGES:
                for direction,pair in [('send','0_1'),('recv','1_0')]:
                    self.write(f'{phase}_{direction}_{pair}_{stage}.txt',[0,17,65535])
            self.write(f'{phase}_send_0_1_selfnodal.txt',[0,17,65535])
    def tearDown(self):self.tmp.cleanup()
    def write(self,name,values):
        (self.root/name).write_text(str(len(values))+'\n'+''.join(str(v)+'\n' for v in values))
    def test_equal(self):self.assertEqual(audit.inspect_packets(self.root)['status'],'PASS')
    def test_inverse_lifting_mismatch(self):
        self.write('restrict_recv_1_0_nodal.txt',[0,18,65535])
        result=audit.inspect_packets(self.root)
        self.assertEqual(result['status'],'MISMATCH')
        self.assertEqual(result['phases']['restrict']['mismatched']['symbols'],0)
        self.assertEqual(result['phases']['restrict']['mismatched']['nodal'],1)
        self.assertEqual(result['mismatches'][0]['first_index'],1)
    def test_predictor_graph_mismatch(self):
        self.write('prolongate_recv_1_0_neighbors.txt',[0,18,65535])
        result=audit.inspect_packets(self.root)
        self.assertEqual(result['status'],'MISMATCH')
        self.assertEqual(result['phases']['prolongate']['mismatched']['symbols'],0)
        self.assertEqual(result['phases']['prolongate']['mismatched']['neighbors'],1)
        self.assertEqual(result['phases']['prolongate']['sender_inverse_mismatched'],0)
    def test_sender_inverse_mismatch(self):
        self.write('prolongate_send_0_1_selfnodal.txt',[0,18,65535])
        result=audit.inspect_packets(self.root)
        self.assertEqual(result['phases']['prolongate']['sender_inverse_mismatched'],1)
    def test_wire_mismatch(self):
        self.write('prolongate_recv_1_0_wire.txt',[1,17,65535])
        self.assertEqual(audit.inspect_packets(self.root)['phases']['prolongate']['mismatched']['wire'],1)
    def test_missing_stage(self):
        (self.root/'restrict_recv_1_0_nodal.txt').unlink()
        self.assertEqual(audit.inspect_packets(self.root)['status'],'INCOMPLETE')
    def test_truncated(self):
        (self.root/'restrict_recv_1_0_nodal.txt').write_text('3\n1\n')
        with self.assertRaises(ValueError):audit.inspect_packets(self.root)
    def test_nonfinite(self):
        self.write('restrict_recv_1_0_values.txt',[float('nan')])
        with self.assertRaises(ValueError):audit.inspect_packets(self.root)

if __name__=='__main__':unittest.main()
