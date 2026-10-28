"""Independent algebraic, source-data, and failure-mode regressions."""
from fractions import Fraction as Q
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from bcov import BCOV, FramedPolynomial, Poly, X, YY_GENUS2_AMBIGUITY, reference_genus2
from bcov.cli import main
from bcov.datasets import iter_amplitudes, load_amplitude
from bcov.patterns import audit_patterns, bernoulli_abs, invert_x, low_x_predictions, pure_x
from bcov.verify_physical_qseries import Series, expected_gv, periods, verify


class PackageTests(unittest.TestCase):
    """Verify public APIs against independently represented source polynomials."""

    def test_source_expansion_and_jets(self) -> None:
        """Atomic jets agree with full expansions in every frame at g=2,3."""
        for g in (2, 3):
            amplitude = load_amplitude(g)
            for frame in ("yy", "hkq", "small_z"):
                polynomial = amplitude.expand(frame).poly.substitute({"u": 0, "v1": 0, "v2": 0, "v3": 0})
                coefficients = pure_x(amplitude, 3*g, frame)
                self.assertEqual(polynomial, sum((c*X**i for i, c in enumerate(coefficients)), Poly()))
                self.assertEqual(pure_x(amplitude, 1, frame), coefficients[:2])
        self.assertEqual(load_amplitude(2).expand("yy").poly, reference_genus2())

    def test_patterns_and_bernoulli(self) -> None:
        """Check independent Bernoulli values and source data through genus 4."""
        self.assertEqual(bernoulli_abs(4), {1: Q(1, 6), 2: Q(1, 30), 3: Q(1, 42), 4: Q(1, 30)})
        for g in (2, 3, 4):
            report = audit_patterns(load_amplitude(g))
            self.assertTrue(report.low_x_matches)
            self.assertTrue(report.denominator_matches)
            self.assertEqual(report.yy_jet[:2], low_x_predictions(g)[g])

    def test_frame_and_inversion_guards(self) -> None:
        """Coordinate inversion is involutive; mixed-frame addition is rejected."""
        vector = (Q(3, 7), Q(-2, 11), Q(5, 13))
        for g in (2, 3, 6):
            self.assertEqual(invert_x(invert_x(vector, g), g), vector)
        p = FramedPolynomial(reference_genus2(), 2, "yy")
        with self.assertRaises(ValueError):
            _ = p + p.to("hkq")
        self.assertEqual(p.to("small_z").to("yy"), p)

    def test_checkpoint_and_unresolved_ambiguity(self) -> None:
        """Unresolved amplitudes fail; exact checkpoints restore the physical seed."""
        solver = BCOV()
        with self.assertRaises(ValueError):
            solver.integrate(2).polynomial()
        solver.set_ambiguity(2, YY_GENUS2_AMBIGUITY)
        solver.integrate(3)
        self.assertEqual(solver.residual(3, solver.amplitudes[3].anomaly), 0)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"checkpoint.json"
            solver.save(path)
            restored = BCOV.load(path)
            self.assertEqual(restored.amplitudes[2].polynomial(), reference_genus2())
            self.assertIsNone(restored.amplitudes[3].ambiguity)

    def test_series_and_mirror_map(self) -> None:
        """Invert exact series and recover the known first mirror-map terms."""
        s = Series(3)
        a = s.trunc([1, -2, 3])
        self.assertEqual(s.mul(a, s.inv(a)), s.c(1))
        self.assertEqual(periods(3)[1][:3], [Q(0), Q(1), Q(-770)])
        for operation in (lambda: s.inv(s.c(0)), lambda: s.exp(s.c(1)), lambda: s.comp(a, s.c(1))):
            with self.assertRaises(ValueError):
                operation()

    def test_missing_data_fails(self) -> None:
        """Incomplete source/GV coverage cannot pass a verification vacuously."""
        with self.assertRaises(ValueError):
            expected_gv(2, 1, [{1: 2875}, {1: 0}, {}])
        with self.assertRaises(ValueError):
            pure_x(load_amplitude(2), -1)
        with self.assertRaises(KeyError):
            load_amplitude(16)
        with tempfile.TemporaryDirectory() as tmp:
            atomic, gv = Path(tmp)/"empty.gz", Path(tmp)/"gv.gz"
            with gzip.open(atomic, "wt") as f:
                f.write("")
            with gzip.open(gv, "wt") as f:
                f.write("genus,degree_1\n0,2875\n1,0\n2,0\n")
            with self.assertRaises(ValueError):
                verify(atomic, gv, 1, 2)
            with gzip.open(atomic, "wt") as f:
                f.write(json.dumps({"genus": 2, "basis": ["wrong frame"], "terms": []}))
            with self.assertRaises(ValueError):
                list(iter_amplitudes(atomic))

    def test_cli_report(self) -> None:
        """The command writes exact, portable rational JSON without package writes."""
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/"report.json"
            self.assertEqual(main(["patterns", "--genus", "2", "--output", str(output)]), 0)
            self.assertEqual(json.loads(output.read_text())[0]["yy_jet"][0], "3125/144")


if __name__ == "__main__":
    unittest.main()
