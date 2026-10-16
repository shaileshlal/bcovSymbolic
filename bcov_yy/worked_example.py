"""A small fully exact BCOV calculation: python worked_example.py."""
from fractions import Fraction as Q

from bcov_yy import BCOV, YY_GENUS2_AMBIGUITY, reference_genus2


def main():
    solver = BCOV()
    print("1. The genus-one one-point seed")
    print("P_1^(1) =", solver.vertex(1,1).expand())

    print("\n2. Form the shared genus-two source, then expand it")
    recipe = solver.source_recipe(2)
    print("Before expansion:", solver.arena.stats())
    source = recipe.expand()
    for power in range(3):
        print(f"[u^{power}] R_2 =", source.coefficient("u",power))

    print("\n3. Recover the three anomaly derivatives")
    for variable, derivative in zip(("v1","v2","v3"),solver.anomaly_gradient(2)):
        print(f"dP_2/d{variable} =", derivative)
    a2 = solver.integrate(2)
    print("Integrated representative:", a2.anomaly)
    print("Ambiguity status:", a2.ambiguity, "(unresolved)")

    print("\n4. Supply published genus-two boundary data")
    p2 = solver.set_ambiguity(2,YY_GENUS2_AMBIGUITY).polynomial()
    assert p2 == reference_genus2()
    assert solver.residual(2) == 0
    print("P_2 matches YY (4.1) exactly; full BCOV residual = 0.")

    print("\n5. Compute the next anomaly polynomial")
    a3 = solver.integrate(3)
    assert a3.anomaly.partial("v3") == Q(-47,3)*p2
    assert solver.residual(3,a3.anomaly) == 0
    print(f"P_3 anomaly: {len(a3.anomaly)} terms, weighted degree {a3.anomaly.degree()}.")
    print("dP_3/dv3 = -47/3 P_2, checked exactly.")
    print("The seven pure-X genus-three coefficients remain unresolved.")
    print("No mirror map or floating point arithmetic was used.")


if __name__ == "__main__":
    main()
