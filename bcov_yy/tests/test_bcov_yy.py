"""Focused exact tests: python -m unittest discover -s tests -v."""
from fractions import Fraction as Q
import json
from pathlib import Path
import random
import tempfile
import unittest

from bcov_yy import (
    Amplitude, BCOV, GENERATORS, NAMES, Poly, RecipeArena, YYRing,
    YY_GENUS2_AMBIGUITY, integrate_gradient, reference_genus2,
    u, v1, v2, v3, X,
)


def sample(seed):
    rng = random.Random(seed)
    terms = {}
    for _ in range(10):
        powers = tuple(rng.randrange(3) for _ in NAMES)
        terms[powers] = Q(rng.randrange(-5, 6), rng.randrange(1, 8))
    return Poly(terms)


class PolynomialTests(unittest.TestCase):
    def test_exactness_cancellation_and_float_rejection(self):
        self.assertEqual(Q(1,3)*v1+Q(2,3)*v1, v1)
        self.assertEqual(len(v1-v1), 0)
        self.assertEqual(Poly.constant(0).degree(), -1)
        for operation in (lambda: Poly.constant(0.1), lambda: v1*0.1,
                          lambda: v1+0.1, lambda: v1/0.1,
                          lambda: v1.evaluate({"v1":0.1})):
            with self.assertRaises(TypeError):
                operation()

    def test_invalid_monomials_and_unsupported_division(self):
        for powers in ((1,0), (-1,0,0,0,0), (1.0,0,0,0,0)):
            with self.assertRaises(ValueError):
                Poly({powers:1})
        with self.assertRaises(ZeroDivisionError):
            _ = v1/0
        with self.assertRaises(ValueError):
            _ = v1**-1
        with self.assertRaises(TypeError):
            _ = v1/(1+X)

    def test_public_terms_are_immutable(self):
        p = v1+v2
        with self.assertRaises(TypeError):
            p.terms[(0,1,0,0,0)] = Q(99)
        self.assertEqual(hash(p), hash(v2+v1))
        self.assertEqual(p, v1+v2)

    def test_arithmetic_against_independent_exact_evaluation(self):
        for seed in range(3):
            p, q, r = sample(seed), sample(seed+10), sample(seed+20)
            values = dict(zip(NAMES, [Q(2,3),Q(-3,5),Q(5,7),Q(7,11),Q(-1,2)]))
            self.assertEqual((p*q).evaluate(values), p.evaluate(values)*q.evaluate(values))
            self.assertEqual(p*(q+r), p*q+p*r)
            self.assertEqual(p*q, q*p)
            self.assertEqual((p+q)-q, p)

    def test_partial_and_primitive(self):
        p, q = sample(1), sample(2)
        for variable in NAMES:
            self.assertEqual(p.integrate(variable).partial(variable), p)
            self.assertEqual((p*q).partial(variable), p.partial(variable)*q+p*q.partial(variable))
        self.assertEqual((v1**3/3).partial("v1"), v1**2)

    def test_simultaneous_substitution_and_coefficients(self):
        self.assertEqual((u+2*v1).substitute({"u":v1,"v1":u}), v1+2*u)
        p = X+u*v1+u**2*v2
        self.assertEqual(p.coefficient("u",0), X)
        self.assertEqual(p.coefficient("u",1), v1)
        self.assertEqual(p.coefficient("u",2), v2)
        self.assertEqual(p.coefficient("u",3), 0)

    def test_weight_bound_and_json_roundtrip(self):
        p = v1**3+v2*X+v3+1
        self.assertEqual(p.degree(), 3)
        self.assertEqual(p.homogeneous_part(3), p-1)
        self.assertEqual(Poly.from_data(json.loads(json.dumps(p.to_data()))), p)


class DifferentialRingTests(unittest.TestCase):
    def setUp(self):
        self.ring = YYRing()

    def test_closed_table_matches_reduced_YY_formulas(self):
        expected = (v2+u*v1-u**2,
                    -v1**2-2*v2-10*X+5*v1*X,
                    -v1*v2+v3,
                    v2**2-24*X-25*v2*X-5*v1*v2*X+10*v3*X,
                    5*X*(X-1))
        self.assertEqual(self.ring.table, expected)
        for generator, derivative in zip(GENERATORS, expected):
            self.assertEqual(self.ring.theta(generator), derivative)
            self.assertLessEqual(derivative.degree(), generator.degree()+1)

    def test_original_AB_Picard_Fuchs_and_special_geometry(self):
        # Check identities in the original variables, not just the transformed table.
        A, B = v1-1-2*u, u
        B2 = v2+u*v1
        B3 = v3-u*v2+5*u*X*(v1-2)
        theta = self.ring.theta
        self.assertEqual(theta(B), B2-B**2)
        self.assertEqual(theta(B2), B3-B*B2)
        B4 = 10*X*B3-35*X*B2+50*X*B-24*X
        self.assertEqual(theta(B3), B4-B*B3)
        A2 = -4*B2-2*A*B-2*B+2*B**2-2*A+10*X*B+5*X*A-5*X-1
        self.assertEqual(theta(A), A2-A**2)

    def test_derivation_and_differential_operator_commutator(self):
        p, q = sample(3), sample(4)
        theta = self.ring.theta
        self.assertEqual(theta(p*q), theta(p)*q+p*theta(q))
        self.assertEqual(theta(p*q)-p*theta(q), theta(p)*q)
        self.assertEqual(theta(Q(7,13)), 0)
        self.assertLessEqual(theta(p).degree(), p.degree()+1)

    def test_z_euler_operator_is_not_partial_z(self):
        self.assertEqual(self.ring.theta_z(X), X*(1-X))
        self.assertEqual(self.ring.theta_z(u*v1), self.ring.theta(u*v1)/-5)

    def test_covariant_derivative_tracks_changing_tensor_rank(self):
        p = v1+X
        first = self.ring.covariant(p,2,0)
        second = self.ring.covariant(first,2,1)
        expected = self.ring.theta(first)+(4*u-v1-5*X)*first
        self.assertEqual(second, expected)
        self.assertNotEqual(second, self.ring.covariant(first,2,0))
        self.assertEqual(self.ring.covariant(Poly.constant(1),0,3), 5*X-3*v1+4*u)


class RecipeTests(unittest.TestCase):
    def test_sharing_and_delayed_expansion(self):
        ring = YYRing()
        arena = RecipeArena(ring)
        a, b = arena.leaf(v1+X), arena.leaf(v2)
        product = a*b
        self.assertEqual(product.node, (b*a).node)
        expression = product+product
        self.assertEqual(arena.stats()["expanded_nodes"], 0)
        self.assertEqual(expression.expand(), 2*(v1+X)*v2)
        self.assertEqual(arena.stats()["operations"]["mul"], 1)
        arena.clear_values()
        self.assertEqual(expression.expand(), 2*(v1+X)*v2)

    def test_iterative_evaluation_avoids_recursion_depth(self):
        arena = RecipeArena(YYRing())
        p = arena.leaf(0)
        for _ in range(1200):
            p = p+1
        self.assertEqual(p.expand(), 1200)

    def test_cross_arena_recipes_rejected(self):
        a, b = RecipeArena(YYRing()), RecipeArena(YYRing())
        with self.assertRaises(ValueError):
            _ = a.leaf(1)+b.leaf(1)


class RecursionTests(unittest.TestCase):
    def test_gradient_reconstruction_and_integrability_failure(self):
        p = v1**2*v2+v1*v3+X*v2**2+X**3
        gradient = [p.partial(v) for v in ("v1","v2","v3")]
        self.assertEqual(integrate_gradient(*gradient), p-X**3)
        with self.assertRaises(ValueError):
            integrate_gradient(v2,Poly(),Poly())

    def test_genus_one_seed_from_original_AB_expression(self):
        solver = BCOV()
        A = v1-1-2*u
        expected = Q(31,3)*(1-u)+Q(5,12)*(X-1)-A/2
        self.assertEqual(solver.vertex(1,1).expand(), expected)

    def test_genus_two_matches_published_polynomial(self):
        solver = BCOV()
        a = solver.set_ambiguity(2, YY_GENUS2_AMBIGUITY)
        self.assertEqual(a.polynomial(), reference_genus2())
        self.assertEqual(solver.residual(2), 0)
        self.assertEqual(a.polynomial().partial("u"), 0)
        self.assertEqual(a.polynomial().degree(), 3)

    def test_all_three_anomaly_derivatives_and_independent_chi_relation(self):
        solver = BCOV()
        p2 = solver.set_ambiguity(2,YY_GENUS2_AMBIGUITY).polynomial()
        a3 = solver.integrate(3)
        for variable,target in zip(("v1","v2","v3"),solver.anomaly_gradient(3)):
            self.assertEqual(a3.anomaly.partial(variable), target)
        self.assertEqual(solver.residual(3,a3.anomaly), 0)
        # From the separate YY partition-function equation (4.7), chi=-200.
        self.assertEqual(a3.anomaly.partial("v3"), Q(-47,3)*p2)
        self.assertEqual(a3.anomaly.coefficient("v3",2), Q(-8225,27))
        self.assertEqual(len(a3.anomaly),55)

    def test_genus_four_nonlinear_chi_relation_and_degree(self):
        solver = BCOV()
        p2 = solver.set_ambiguity(2, YY_GENUS2_AMBIGUITY).polynomial()
        # Explicit formal representative, not claimed to be physical genus-three data.
        p3 = solver.set_ambiguity(3,()).polynomial()
        a4 = solver.integrate(4)
        self.assertEqual(a4.anomaly.partial("v3"), Q(-82,3)*p3+2*p2**2)
        self.assertEqual(a4.anomaly.degree(),9)
        self.assertEqual(a4.anomaly.partial("u"),0)
        self.assertEqual(solver.residual(4,a4.anomaly),0)

    def test_unresolved_ambiguity_blocks_continuation(self):
        solver = BCOV()
        a2 = solver.integrate(2)
        with self.assertRaisesRegex(ValueError,"unresolved"):
            a2.polynomial()
        with self.assertRaisesRegex(ValueError,"unresolved"):
            solver.integrate(3)
        solver.set_ambiguity(2,())
        self.assertEqual(solver.integrate(3).anomaly.degree(),6)

    def test_lower_ambiguity_propagates_and_stale_vertices_are_invalidated(self):
        a, b = BCOV(), BCOV()
        a.set_ambiguity(2,())
        old = a.source_recipe(3).expand()
        a.set_ambiguity(2,(Q(7,11),))
        b.set_ambiguity(2,(Q(7,11),))
        self.assertNotEqual(a.source_recipe(3).expand(),old)
        self.assertEqual(a.integrate(3).anomaly,b.integrate(3).anomaly)
        with self.assertRaisesRegex(ValueError,"invalidate"):
            a.set_ambiguity(2,())

    def test_holomorphic_kernel_and_degree_guard(self):
        solver = BCOV()
        a = solver.integrate(2)
        source = solver.source_recipe(2).expand()
        p = a.with_ambiguity((Q(1,7),Q(2,9),0,Q(3,11))).polynomial()
        self.assertEqual(solver.residual(2,p),0)
        self.assertEqual(source,solver.source_recipe(2).expand())
        with self.assertRaises(ValueError):
            a.with_ambiguity((0,0,0,0,1))
        with self.assertRaises(ValueError):
            Amplitude(2,u)
        with self.assertRaises(ValueError):
            Amplitude(2,v1**4)
        with self.assertRaises(ValueError):
            Amplitude(2,X)

    def test_checkpoint_roundtrip_resume_and_tamper_detection(self):
        solver = BCOV()
        solver.set_ambiguity(2,YY_GENUS2_AMBIGUITY)
        solver.integrate(3)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"checkpoint.json"
            solver.save(path)
            restored = BCOV.load(path)
            self.assertEqual(restored.amplitudes[2].polynomial(),reference_genus2())
            self.assertEqual(restored.amplitudes[3].anomaly,solver.amplitudes[3].anomaly)
            self.assertIsNone(restored.amplitudes[3].ambiguity)
            restored.set_ambiguity(3,())
            self.assertEqual(restored.integrate(4).anomaly.degree(),9)
            data = json.loads(path.read_text())
            data["amplitudes"][1]["anomaly"][0][1] = "9999999999"
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,"anomaly equation"):
                BCOV.load(path)

    def test_compaction_preserves_results(self):
        solver = BCOV()
        solver.set_ambiguity(2,YY_GENUS2_AMBIGUITY)
        a3 = solver.integrate(3)
        solver.compact_recipes()
        self.assertEqual(solver.arena.stats()["nodes"],0)
        self.assertEqual(solver.residual(3,a3.anomaly),0)


if __name__ == "__main__":
    unittest.main()
