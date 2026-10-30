import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_emi_bddc_convergence import read_history


class HistoryChecks(unittest.TestCase):
    def read(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "history.csv"
            with path.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("time_step", "iteration", "step_length", "residual_norm",
                                 "rhs_norm", "relative_residual", "target"))
                writer.writerows(rows)
            return read_history(path)

    def test_initial_and_updated_rows(self):
        rows = ((1, 0, "", 2, 2, 1, 1e-7), (1, 1, 0.2, 0.5, 2, 0.25, 1e-7))
        history = self.read(rows)
        self.assertIsNone(history[0].step_length)
        self.assertEqual(history[1].relative_residual, 0.25)

    def test_zero_rhs_absolute_fallback(self):
        self.assertEqual(self.read(((2, 0, "", 0, 0, 0, 1e-7),))[0].rhs_norm, 0)

    def test_inconsistent_or_invalid_data_rejected(self):
        bad_rows = (
            ((1, 0, ""),),
            ((1, 0, 0.2, 1, 1, 1, 1e-7),),
            ((1, 0, "", 1, 2, 1, 1e-7),),
            ((1, 0, "", "nan", 1, "nan", 1e-7),),
            ((1, 0, "", 1, 1, 1, 0),),
            ((1, 0, "", 1, 1, 1, 1e-7), (1, 2, 0.2, 0.5, 1, 0.5, 1e-7)),
            ((1, 0, "", 1, 1, 1, 1e-7), (1, 1, "", 0.5, 1, 0.5, 1e-7)),
            ((1, 0, "", 1, 1, 1, 1e-7), (1, 1, 0.2, 1, 2, 0.5, 1e-7)),
        )
        for rows in bad_rows:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                self.read(rows)


if __name__ == "__main__":
    unittest.main()
