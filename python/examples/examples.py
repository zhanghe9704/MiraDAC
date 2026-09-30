"""Port of examples/examples.cc: NDA arithmetic, substitution and composition."""
import miradac as da


def main() -> None:
    print("Print out the base vector.\n")
    for i in range(3):
        print(da.base[i], end="")

    print("Fundamental calculations of DA vectors.\n")
    x = 1.0 + da.base[0] + 2.0 * da.base[1] + 5.0 * da.base[2]
    print(x, end="")

    y = da.exp(x)
    print(y, end="")

    print("Substitute a number for a base.\n")
    z = da.NDA()
    da.substitute(y, 0, 1.0, z)
    print(z, end="")

    print("Substitute a DA vector for a base.\n")
    da.substitute(y, 0, x, z)
    print(z, end="")

    print("Substitute multiple DA vectors for bases at once.\n")
    lv = da.NDAList([da.sin(x), da.cos(x)])
    idx = [0, 1]
    da.substitute(y, idx, lv, z)
    print(z, end="")

    print("The norm of z is", z.norm())
    print("The weighted norm of z is", z.weighted_norm(0.1))

    print("Bunch processing for substitutions.\n")
    lx = da.NDAList([x, y, da.sinh(x)])
    ly = da.NDAList([da.NDA() for _ in range(3)])
    da.substitute(lx, idx, lv, ly)
    for v in ly:
        print(v, end="")

    print("Composition of DA vectors with numbers.\n")
    ln = da.compose(lx, [0.1, 2, 1])
    print(*ln, "\n")

    print("Composition of DA vectors with DA vectors.\n")
    lu = da.NDAList([da.sin(x), da.cos(x), da.tan(x)])
    da.compose(lx, lu, ly)
    for v in ly:
        print(v, end="")

    # Output a DA vector to file.
    with open("da_output.txt", "w") as f:
        f.write(str(z))

    z /= 1e5
    print("z:")
    print(z, end="")
    print("norm of z:", z.norm())

    print("DA eps:", da.get_eps())
    da.set_eps(1e-20)
    print("Reset DA eps:", da.get_eps())


if __name__ == "__main__":
    da.init(order=4, nvars=3, pool_size=100)
    main()  # every DA vector is freed when main() returns, before clear()
    da.clear()
