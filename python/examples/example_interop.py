"""Port of examples/example_interop.cc: NDA * SDA -> SDA, then evaluate -> NDA."""
import sys

import miradac as da


def main() -> int:
    print("=== MiraDAC Interop Example ===")
    print(f"order={da.current_order()}  nv={da.nvars()}  full_len={da.full_length()}\n")

    # 1. Numerical DA vector: nx = x (the first base variable).
    nx = da.var(0)
    print("NDA: nx = var(0)")
    print(nx)

    # 2. Symbolic DA vector: sy = alpha * y.
    alpha = da.Expr("alpha")
    sy = alpha * da.var(1)
    print("\nSDA: sy = alpha * var(1)")
    print(sy)

    # 3. Mixed multiplication: NDA * SDA -> SDA (the NDA is promoted).
    result_sda = nx * sy
    print("\nSDA result = nx * sy  (NDA * SDA -> SDA, via promote):")
    print(result_sda)

    # 4. Evaluate: substitute alpha = 3.0, get an NDA.
    result_nda = da.evaluate(result_sda, {alpha: 3.0})
    print("\nNDA result = evaluate(result_sda, {alpha: 3.0}):")
    print(result_nda)

    # The coefficient of x*y should be 3.0.
    coeff = result_nda.element([1, 1])
    print(f"\nCoefficient of x*y = {coeff:.6f}  (expected 3.0)")
    if abs(coeff - 3.0) < 1e-12:
        print("PASSED: coefficient matches 3.0")
        return 0
    print("FAILED: coefficient mismatch!")
    return 1


if __name__ == "__main__":
    da.init(order=4, nvars=2, pool_size=500, table=True)
    status = main()
    da.clear()
    sys.exit(status)
