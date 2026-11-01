import csv
import json
from pathlib import Path
import tempfile
import unittest
from run_emi_precision_followup import can_extend, CONTROLS, MODE

class FollowupTest(unittest.TestCase):
    def test_cap_extension_requires_valid_decreasing_run_and_controls(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);folder=root/"10cells"/MODE;folder.mkdir(parents=True)
            records={m:dict(status="PASS") for m in CONTROLS}
            records[MODE]=dict(status="FAIL",saved_trace_prefix_rows=301)
            (folder/"run.log").write_text("Step 1: iteration limit")
            def save(vals):
                (root/"summary.json").write_text(json.dumps(records))
                with (folder/"bddc-fgmres-step1.csv").open("w",newline="") as f:
                    w=csv.writer(f);w.writerow(["iteration","true_relative_residual"])
                    w.writerows(enumerate(vals))
            save([0.99**i for i in range(501)])
            self.assertTrue(can_extend(root,500)[0])
            save([0.99999**i for i in range(501)])
            self.assertTrue(can_extend(root,500)[0])
            save([1.0]*501);self.assertFalse(can_extend(root,500)[0])
            save([0.99**i for i in range(500)]);self.assertFalse(can_extend(root,500)[0])
            records["full-fgmres16"]["status"]="FAIL"
            save([0.99**i for i in range(501)]);self.assertFalse(can_extend(root,500)[0])
            records["full-fgmres16"]["status"]="PASS";records[MODE]["status"]="PASS"
            save([0.99**i for i in range(501)]);self.assertFalse(can_extend(root,500)[0])
            records[MODE]["status"]="FAIL";records[MODE].pop("saved_trace_prefix_rows")
            save([0.99**i for i in range(501)]);self.assertFalse(can_extend(root,500)[0])

if __name__=="__main__":unittest.main()
