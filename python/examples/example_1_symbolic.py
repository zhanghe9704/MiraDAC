"""Port of examples/example_1_symbolic.cc: SDA fundamental operations."""
import miradac as da


def main() -> None:
    x, y, z = da.symbols("x y z")
    b = da.base

    # An int stays an exact integer coefficient, a float becomes a float.
    da1 = 1 + (1 + x) * b[0] + y * b[1] + (z - 0.5) * b[2]
    print("Symbolic DA vector 1:", da1, sep="\n", end="")

    da2 = 3.3 + (0.5 + x) * b[0] + y * y * b[1] + (x + z + 1.1) * b[2]
    print("Symbolic DA vector 2:", da2, sep="\n", end="")

    print("Summation: vec 1 + vec 2:")
    print(da1 + da2, end="")

    print("Multiplication: vec 1 * vec 2:")
    print(da1 * da2, end="")

    print("Inverse of vec 1:")
    print(1 / da1, end="")

    print("Square root of vec 1:")
    print(da.sqrt(da1), end="")


if __name__ == "__main__":
    da.init(order=4, nvars=3, pool_size=400, table=True)
    main()
    da.clear()
