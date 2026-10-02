# Port of examples/examples.cc: NDA arithmetic, substitution and composition.
using MiraDAC
using LinearAlgebra: norm

function main()
    println("Print out the base vector.\n")
    foreach(i -> display(davar(i)), 1:3)

    println("Fundamental calculations of DA vectors.\n")
    x = 1.0 + davar(1) + 2.0 * davar(2) + 5.0 * davar(3)
    display(x)

    y = exp(x)
    display(y)

    println("Substitute a number for a base.\n")
    z = substitute(y, 1, 1.0)
    display(z)

    println("Substitute a DA vector for a base.\n")
    z = substitute(y, 1, x)
    display(z)

    println("Substitute multiple DA vectors for bases at once.\n")
    lv = NDAList([sin(x), cos(x)])
    idx = [1, 2]
    z = substitute(y, idx, lv)
    display(z)

    println("The norm of z is ", norm(z))
    # C++ also prints z.weighted_norm(0.1); the Julia API (plan A.6) has no weighted norm.

    println("Bunch processing for substitutions.\n")
    lx = NDAList([x, y, sinh(x)])
    ly = substitute(lx, idx, lv)
    foreach(display, ly)

    println("Composition of DA vectors with numbers.\n")
    println(join(compose(lx, [0.1, 2, 1]), ' '), "\n")

    println("Composition of DA vectors with DA vectors.\n")
    lu = NDAList([sin(x), cos(x), tan(x)])
    ly = compose(lx, lu)
    foreach(display, ly)

    # Output a DA vector to file.
    write("da_output.txt", repr(MIME"text/plain"(), z))

    z /= 1e5
    println("z:")
    display(z)
    println("norm of z: ", norm(z))

    println("DA eps: ", get_eps())
    set_eps!(1e-20)
    println("Reset DA eps: ", get_eps())
end

init!(4, 3, 100)
main()
clear!()
