import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from collect_emi_scaling_smoke import collect


def log(steps=1, ranks=8, status="converged", rel="7e-8", compression=0):
    rows="\n".join(f"BDDC global residual step {i}: abs=1e-10 rel={rel} "
                   f"target=1e-07 status={status} iterations=32" for i in range(1,steps+1))
    return f"""cells: 9793
dofs: 3164
BDDC subdomains: 16, shared dofs: 359
BDDC coarse interfaces: corners=1, edges=42, faces=32
{rows}
Compression report (BDDC):
  raw sent bytes: 10000
  transmitted sent bytes: {4000 if compression else 10000}
Profile (wall seconds, EMI end-to-end phases; {ranks} rank(s)):
  solver wrapper (including setup)              1.0 1.1 1.2
total cpu-time:  2.300000s wall, 3s user + 1s system = 4s CPU
EMI model completed
"""


class Metrics(unittest.TestCase):
    def test_valid(self):
        result=collect(log(steps=2,compression=1),2,8,0)
        self.assertEqual(result["status"],"converged")
        self.assertEqual(result["iterations_per_step"],[32,32])
        self.assertEqual(result["reported_transfer_savings_percent"],60)
        self.assertEqual(result["solver_wrapper_max_seconds"],1.2)

    def test_explicit_tolerance(self):
        text = log(rel="7e-7").replace("target=1e-07", "target=1e-06")
        result = collect(text, 1, 8, 0, expected_tolerance=1e-6)
        self.assertEqual(result["status"], "converged")
        self.assertEqual(result["requested_relative_tolerance"], 1e-6)
        self.assertEqual(collect(text, 1, 8, 0)["status"], "failed")
        for bad in (text.replace("rel=7e-7", "rel=2e-6"),
                    text.replace("target=1e-06", "target=1e-05")):
            self.assertEqual(collect(bad, 1, 8, 0, 1e-6)["status"], "failed")
        for invalid in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                collect(text, 1, 8, 0, invalid)

    def test_reject_incomplete_or_invalid(self):
        for text,steps,ranks,rc in (
            (log(status="not-converged"),1,8,0),
            (log(rel="2e-7"),1,8,0),
            (log(rel="nan"),1,8,0),
            (log(rel="1e999"),1,8,0),
            (log(),2,8,0),
            (log(),1,4,0),
            (log(),1,8,137),
            (log().replace("EMI model completed",""),1,8,0),
            (log().replace("target=1e-07","target=1e-05"),1,8,0),
            (log()+"  raw sent bytes: 10000\n",1,8,0),
        ):
            with self.subTest(text=text,steps=steps,ranks=ranks,rc=rc):
                result=collect(text,steps,ranks,rc)
                self.assertEqual(result["status"],"failed")
                json.dumps(result,allow_nan=False)

    def test_expansion_is_a_valid_measurement(self):
        result=collect(log(compression=1).replace("bytes: 4000","bytes: 11000"),1,8,0)
        self.assertEqual(result["status"],"converged")
        self.assertLess(result["reported_transfer_savings_percent"],0)


class Runner(unittest.TestCase):
    def test_flags_isolation_and_failure(self):
        source=Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory(prefix="emi-scaling-control-") as folder:
            root=Path(folder)
            repo=root/"repo"
            (repo/"tests").mkdir(parents=True)
            for name in ("run_emi_scaling_smoke.sh","collect_emi_scaling_smoke.py"):
                shutil.copyfile(source/name,repo/"tests"/name)
            inputs=root/"meshes"
            inputs.mkdir()
            for case in ("weak_8_h16","strong_h16"):
                for suffix in ("_mesh.vtu","_list_extracellular.txt",
                               "_list_intracellular.txt","_early_excited.txt"):
                    (inputs/(case+"_kaskade"+suffix)).touch()
            model=repo/"emiModel"
            model.write_text("#!/bin/sh\nexit 0\n")
            model.chmod(0o755)
            fake=root/"mpi"
            fake.write_text("#!/usr/bin/env python3\n"+'''
import json,os,sys
from pathlib import Path
if '--version' in sys.argv:
    print('fake MPI'); sys.exit(0)
args=sys.argv[1:]
def value(key): return args[args.index(key)+1]
with open(os.environ['CAPTURE'],'a') as f: f.write(json.dumps(args)+'\\n')
steps=int(value('--maximumNumberOfTimeSteps'))
ranks=int(value('-np'))
comp=int(value('--bddcCompression'))
run=Path(value('--dir'))
bad=comp and os.environ.get('FAIL_FULL')=='1'
print('cells: 9793\\ndofs: 3164')
print('BDDC subdomains: 16, shared dofs: 359')
print('BDDC coarse interfaces: corners=1, edges=42, faces=32')
for i in range(1,steps+1):
    print(f"BDDC global residual step {i}: abs=1e-10 rel={2e-7 if bad else 7e-8} target=1e-07 status={'not-converged' if bad else 'converged'} iterations=32")
print('Compression report (BDDC):\\n  raw sent bytes: 10000')
print('  transmitted sent bytes: '+str(4000 if comp else 10000))
print(f'Profile (wall seconds, EMI end-to-end phases; {ranks} rank(s)):')
print('  solver wrapper (including setup) 1.0 1.1 1.2')
print('total cpu-time:  2.300000s wall, 3s user + 1s system = 4s CPU')
print('EMI model completed')
(run/('emiBDDCLastCompression.vtu' if comp else 'emiBDDCLast.vtu')).touch()
sys.exit(17 if bad else 0)
''')
            fake.chmod(0o755)
            (repo/"compare_emi_vtu.py").write_text(
                "from pathlib import Path\nimport sys\n"
                "assert len(sys.argv) == 3\n"
                "assert all(Path(p).is_file() for p in sys.argv[1:])\n"
                "assert Path(sys.argv[1]).name == 'emiBDDCLastCompression.vtu'\n"
                "assert Path(sys.argv[2]).name == 'emiBDDCLast.vtu'\n"
                "print('relative L2 error: 2e-9')\nprint('max absolute error: 1e-8')\n")
            capture=root/"args.jsonl"
            # The production container has no hostname utility. A deliberately
            # failing shadow catches accidental reintroduction of that dependency.
            fake_bin=root/"bin"
            fake_bin.mkdir()
            hostname=fake_bin/"hostname"
            hostname.write_text("#!/bin/sh\necho HOSTNAME_PROBE >&2\nexit 127\n")
            hostname.chmod(0o755)
            env=dict(os.environ,MPIEXEC=str(fake),CAPTURE=str(capture),MESH_DIR=str(inputs),
                     OUT_ROOT=str(root/"success"),STEPS="1",RANKS="8",BITS="12",
                     CASES="weak_8_h16 strong_h16",MAX_ITERATIONS="300",
                     PATH=str(fake_bin)+os.pathsep+os.environ["PATH"])
            completed=subprocess.run(["bash",str(repo/"tests/run_emi_scaling_smoke.sh")],
                                     env=env,capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stderr)
            metadata=(root/"success/environment.txt").read_text()
            self.assertIn("hostname=",metadata)
            self.assertNotIn("HOSTNAME_PROBE",metadata)
            args=[json.loads(line) for line in capture.read_text().splitlines()]
            self.assertEqual(len(args),4)
            for row in args:
                def value(key): return row[row.index(key)+1]
                for flag in ("--bddcGraphLifting","--bddcGraphLiftingAfterQuantization",
                             "--bddcZigzag","--bddcBitlength","--bddcHuffman"):
                    self.assertEqual(value(flag),"1")
                self.assertEqual(value("--bddcUseCg"),"1")
                self.assertEqual(value("--bddcRelativeTolerance"),"1e-7")
                self.assertEqual(value("--bddcConditionEstimate"),"0")
                self.assertEqual(value("--vtk"),"1")
            failure_env=dict(env,OUT_ROOT=str(root/"failure"),FAIL_FULL="1")
            failed=subprocess.run(["bash",str(repo/"tests/run_emi_scaling_smoke.sh")],
                                  env=failure_env,capture_output=True,text=True)
            self.assertNotEqual(failed.returncode,0)
            metrics=json.loads((root/"failure/weak_8_h16/full-bits12/metrics.json").read_text())
            self.assertEqual(metrics["status"],"failed")
            self.assertFalse((root/"failure/weak_8_h16/comparison.log").exists())
            self.assertFalse((root/"failure/strong_h16").exists())


if __name__ == "__main__":
    unittest.main()
