"""Small reproducible formal scaling check, not a physical high-genus solution."""
import argparse
import json
from time import perf_counter

from bcov_yy import BCOV, YY_GENUS2_AMBIGUITY


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", type=int, default=8)
    args = parser.parse_args()
    if args.genus < 2:
        parser.error("--genus must be >= 2")
    solver = BCOV()
    print("# Genus 2: reference boundary data. Later genera: explicit zero ambiguity.")
    for genus in range(2,args.genus+1):
        start = perf_counter()
        a = solver.integrate(genus)
        if solver.residual(genus,a.anomaly) != 0:
            raise AssertionError("Nonzero BCOV residual")
        elapsed = perf_counter()-start
        coefficient_bits = max(max(abs(c.numerator).bit_length(),c.denominator.bit_length())
                               for c in a.anomaly.terms.values())
        print(json.dumps({"genus":genus,"anomaly_terms":len(a.anomaly),
                          "weighted_degree":a.anomaly.degree(),"seconds":round(elapsed,6),
                          "max_coefficient_bits":coefficient_bits,"residual_zero":True}))
        if genus < args.genus:
            solver.set_ambiguity(genus,YY_GENUS2_AMBIGUITY if genus == 2 else ())
        # Release temporary expansions between genera to bound cache overhead.
        solver.compact_recipes()


if __name__ == "__main__":
    main()
