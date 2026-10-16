from fractions import Fraction as Q
from pathlib import Path
import tempfile
import unittest

from bcov_yy import BCOV, YY_GENUS2_AMBIGUITY, YYRing, X, v1, v2, set_resource_limits, ResourceLimitError
from exact_series import Jet, evaluate_poly
from period_frames import frame, conifold_periods, conifold_action
from physical_boundaries import (solve_boundary, seed_gv, record_gv, required_gv,
                                constant_map, conifold_coefficient, ExactLinearSystem,
                                gv_kernel, q_coefficients, gw_from_gv)
from monodromy import validate
from run_physical import plan, run


class SeriesAndFrameTests(unittest.TestCase):
    def test_precision_and_laurent_inversion(self):
        t=Jet({1:1},12)
        p=(t*(1-t)).inverse()
        self.assertEqual(p.valuation,-1)
        self.assertEqual(p.precision,10)
        self.assertEqual([p.coeff(i) for i in range(-1,6)],[1]*7)
        with self.assertRaises(ValueError): p.coeff(10)
        self.assertEqual((p*t*(1-t)-1).c,{})

    def test_local_monodromies(self):
        self.assertTrue(validate())

    def test_monodromy_image_period(self):
        f,nu=conifold_periods(14)
        self.assertEqual([nu.coeff(i) for i in range(1,4)],[1,Q(7,10),Q(41,75)])
        # Independent PF action on the generated analytic series.
        for k in range(10):
            c=sum((a*conifold_action(Q(j)).get(k-j,(0,0))[0] for j,a in nu.c.items()),Q(0))
            self.assertEqual(c,0)
        # The obstruction to an analytic completion of nu log(delta) vanishes.
        c=sum((a*conifold_action(Q(j)).get(-1-j,(0,0))[1] for j,a in nu.c.items()),Q(0))
        self.assertEqual(c,0)

    def test_ring_closure_in_all_frames(self):
        ring=YYRing()
        n=32
        for name in ['mum','conifold','orbifold']:
            values,_,_=frame(name,n)
            theta=(lambda a:-5*a.euler()) if name=='mum' else (
                (lambda a:5*(1-Jet({1:1},n))*a.derivative()) if name=='conifold' else (lambda a:5*a.euler()))
            for value,relation in zip(values,ring.table):
                residual=theta(value)-evaluate_poly(relation,values)
                for k in range(residual.valuation,min(10,residual.precision)):
                    self.assertEqual(residual.coeff(k),0,(name,k))

    def test_q_lagrange_inversion(self):
        _,_,q=frame('mum',20)
        self.assertEqual(q_coefficients(q+3*q*q,5),[0,1,3,0,0,0])

    def test_gv_seed_and_kernel(self):
        gv=seed_gv(4)
        self.assertEqual(gv[0,1],2875)
        self.assertEqual(gv[0,4],242467530000)
        self.assertEqual(gv[1,3],609250)
        self.assertEqual(gv_kernel(2,0),Q(1,240))
        self.assertEqual(gv_kernel(3,0),Q(1,6048))

    def test_constant_and_gap_signs(self):
        self.assertEqual(constant_map(2),Q(-5,144))
        self.assertEqual(constant_map(3),Q(5,36288))
        self.assertEqual(conifold_coefficient(2),Q(1,1200))
        self.assertEqual(conifold_coefficient(3),Q(1,25200))


class PhysicalBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s=BCOV()
        cls.gv=seed_gv(4)
        cls.diags={}
        for g in range(2,5):
            a=cls.s.integrate(g)
            vector,diag=solve_boundary(g,a.anomaly,cls.gv)
            cls.s.set_ambiguity(g,vector)
            record_gv(g,cls.s.amplitudes[g].polynomial(),4,cls.gv)
            cls.diags[g]=diag

    def test_boundary_derives_reference_genus_two(self):
        self.assertEqual(self.s.amplitudes[2].ambiguity,YY_GENUS2_AMBIGUITY)

    def test_higher_genus_integer_gv(self):
        self.assertEqual(self.gv[2,4],534750)
        self.assertEqual(self.gv[3,4],8625)
        self.assertEqual(self.gv[4,4],0)
        self.assertTrue(all(c.denominator==1 for c in self.gv.values()))

    def test_physical_anomaly_residuals(self):
        for g in range(2,5): self.assertEqual(self.s.residual(g),0)

    def test_jet_precision_stability(self):
        for g in [2,3]:
            v,_=solve_boundary(g,self.s.amplitudes[g].anomaly,self.gv,precision=50)
            self.assertEqual(v,self.s.amplitudes[g].ambiguity)

    def test_full_boundary_coefficients(self):
        for g in range(2,5):
            poly=self.s.amplitudes[g].polynomial()
            values,k,t=frame('conifold',50)
            gap=evaluate_poly(poly,values)*k**(g-1)-conifold_coefficient(g)*t**(2-2*g)
            self.assertTrue(all(gap.coeff(j)==0 for j in range(2-2*g,0)))
            values,k,_=frame('orbifold',30)
            orb=evaluate_poly(poly,values)*k**(g-1)
            self.assertTrue(all(orb.coeff(j)==0 for j in range((3*g+1)//5)))

    def test_rank_defect_and_inconsistency(self):
        s=ExactLinearSystem(2)
        s.add([1,1],2,'first')
        with self.assertRaises(ValueError): s.solve()
        with self.assertRaises(ValueError): s.add([1,1],3,'contradiction')

    def test_genus_51_exception_and_target_53(self):
        self.assertEqual(required_gv(51,20),175)
        self.assertEqual(required_gv(52,20),0)
        self.assertEqual(required_gv(53,20),0)
        self.assertIsNone(required_gv(50,20))
        pp=plan()
        self.assertEqual(pp[-1]['genus'],53)
        self.assertEqual(pp[-1]['unknowns'],157)
        self.assertEqual(pp[-1]['gv_max_degree'],20)
        for p in pp:
            self.assertEqual(p['conifold']+p['orbifold']+p['mum'],p['unknowns'])

    def test_resource_limit(self):
        try:
            set_resource_limits(0,1)
            with self.assertRaises(ResourceLimitError): (1+X)*(v1+v2)
        finally: set_resource_limits()

    def test_resource_stop_preserves_physical_checkpoint(self):
        import json
        with tempfile.TemporaryDirectory() as d:
            try:
                with self.assertRaises(ResourceLimitError):
                    run(3,d,max_products=20)
                status=json.loads((Path(d)/'run.json').read_text())
                self.assertEqual(status['status'],'stopped')
                self.assertEqual(status['completed_genus'],2)
                set_resource_limits()
                cp=BCOV.load(Path(d)/'amplitudes.json')
                self.assertEqual(cp.amplitudes[2].ambiguity,YY_GENUS2_AMBIGUITY)
            finally: set_resource_limits()

    def test_resume_rejects_formal_ambiguity(self):
        with tempfile.TemporaryDirectory() as d:
            s=BCOV()
            s.integrate(2)
            s.set_ambiguity(2,())
            s.save(Path(d)/'amplitudes.json')
            with self.assertRaisesRegex(ValueError,'physical boundaries'):
                run(3,d)
        set_resource_limits()

    def test_symbolic_mu_resume(self):
        with tempfile.TemporaryDirectory() as d:
            s,_=run(2,d,mu='mu')
            t,_=run(3,d,mu='mu')
            self.assertEqual(s.amplitudes[2].polynomial(),t.amplitudes[2].polynomial())
            self.assertEqual(t.amplitudes[3].polynomial(),self.s.amplitudes[3].polynomial())
            self.assertIn('"mu": "mu"',(Path(d)/'run.json').read_text())
        set_resource_limits()


if __name__=='__main__': unittest.main()
