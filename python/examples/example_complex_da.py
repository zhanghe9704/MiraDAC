"""Port of examples/example_complex_da.cc: CNDA arithmetic and composition."""
import math

import miradac as da


def main() -> None:
    x1 = da.var(0) + 2.0 * da.var(1) + 3.0 * da.var(2)
    x2 = da.sin(x1)
    x1 = da.cos(x1)

    x3 = 0.5 * da.var(0) + 4.0 * da.var(1) + 2.7 * da.var(2)
    x4 = da.sin(x3)
    x3 = da.cos(x3)

    # CNDA(re, im) is re + i*im; there is no NDA + CNDA operator.
    y1 = da.CNDA(x1, x2)
    y2 = da.CNDA(x3, x4)

    print("y1:", y1, sep="\n")
    print("y2:", y2, sep="\n")

    print("y1+y2:", y1 + y2, sep="\n")
    print("y1-y2:", y1 - y2, sep="\n")
    print("y1*y2:", y1 * y2, sep="\n")
    print("y1/y2:", y1 / y2, sep="\n")

    mmap = da.NDAList([x1, x2])
    nmap = [4.2 + 0.3j, 1 / 3.0 + math.sqrt(2.0) * 1j, math.sin(0.7) + math.cos(0.4) * 1j]

    omap = da.compose(mmap, nmap)
    print("Composition of DA vectors with complex numbers.\n")
    for o in omap:
        print(o)

    cnmap = da.CNDAList([y1, y2, y1 * y2])
    comap = da.CNDAList([da.CNDA(), da.CNDA()])
    da.cd_composition(mmap, cnmap, comap)
    print("Composition of DA vectors with complex DA vectors.\n")
    for c in comap:
        print(c)

    cmmap = da.CNDAList([da.CNDA(x1, da.exp(x1)), da.CNDA(x2, da.exp(x2))])
    da.cd_composition(cmmap, cnmap, comap)
    print("Composition of complex DA vectors with complex DA vectors.\n")
    for c in comap:
        print(c)

    mmap.append(x1 + 0.33 * x2)
    da.cd_composition(cmmap, mmap, comap)
    print("Composition of complex DA vectors with DA vectors.\n")
    for c in comap:
        print(c)


if __name__ == "__main__":
    da.init(order=4, nvars=3, pool_size=400)
    main()
    da.clear()
