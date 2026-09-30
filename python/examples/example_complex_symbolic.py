"""Port of examples/example_complex_symbolic.cc: symbolic complex DA (CSDA).

CSDA keeps the imaginary unit outside the coefficients, so the real and
imaginary parts stay separable and evaluate() returns a CNDA.
"""
import math
import sys

import miradac as da


def main() -> int:
    a, b = da.symbols("a b")

    # z = (1.5 + a*x1 + 0.1*x2) + i*(0.5 + b*x1 - 0.05*x2). The constant parts
    # are numeric so that exp, sqrt, ... can branch on them.
    re = 1.5 + a * da.svar(0) + 0.1 * da.svar(1)
    im = 0.5 + b * da.svar(0) - 0.05 * da.svar(1)
    z = da.CSDA(re, im)

    print("=== Symbolic Complex DA example ===\n")
    print("z = (1.5 + a*x1 + 0.1*x2) + i*(0.5 + b*x1 - 0.05*x2)\n")

    z2 = z * z
    print("z^2 real constant part:", float(z2.real.con))  # 1.5^2 - 0.5^2 = 2.0
    print("z^2 imag constant part:", float(z2.imag.con))  # 2*1.5*0.5 = 1.5

    ez = da.exp(z)
    print(f"\nexp(z) real constant part: {float(ez.real.con)} "
          f"(expected ~{math.exp(1.5) * math.cos(0.5)})")
    print(f"exp(z) imag constant part: {float(ez.imag.con)} "
          f"(expected ~{math.exp(1.5) * math.sin(0.5)})")

    re_ez, im_ez = ez.real, ez.imag
    print("\nReal part of exp(z) at order 0:", float(re_ez.con))
    print("Imag part of exp(z) at order 0:", float(im_ez.con))

    # Evaluate at a = 0.3, b = 0.2 and compare with the numeric computation.
    vals = {a: 0.3, b: 0.2}
    ez_num = da.evaluate(ez, vals)

    n = da.CNDA(1.5 + 0.3 * da.var(0) + 0.1 * da.var(1),
                0.5 + 0.2 * da.var(0) - 0.05 * da.var(1))
    ok = da.compare_cd_vectors(ez_num, da.exp(n), 1e-12)
    print("\nEvaluation cross-check (exp):", "PASS" if ok else "FAIL")

    sz_num = da.evaluate(da.sqrt(z), vals)
    ok2 = da.compare_cd_vectors(sz_num, da.sqrt(n), 1e-12)
    print("Evaluation cross-check (sqrt):", "PASS" if ok2 else "FAIL")

    if not (ok and ok2):
        return 1
    print("\nAll CSDA operations completed successfully.")
    return 0


if __name__ == "__main__":
    da.init(order=3, nvars=2, pool_size=300)
    status = main()
    da.clear()
    sys.exit(status)
