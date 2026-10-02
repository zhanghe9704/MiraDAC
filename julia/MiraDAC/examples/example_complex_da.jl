# Port of examples/example_complex_da.cc: CNDA arithmetic and composition.
using MiraDAC

function main()
    x1 = davar(1) + 2.0 * davar(2) + 3.0 * davar(3)
    x2 = sin(x1)
    x1 = cos(x1)

    x3 = 0.5 * davar(1) + 4.0 * davar(2) + 2.7 * davar(3)
    x4 = sin(x3)
    x3 = cos(x3)

    # CNDA(re, im) is re + i*im; there is no NDA + CNDA operator.
    y1 = CNDA(x1, x2)
    y2 = CNDA(x3, x4)

    for (name, v) in ("y1" => y1, "y2" => y2, "y1+y2" => y1 + y2, "y1-y2" => y1 - y2,
                      "y1*y2" => y1 * y2, "y1/y2" => y1 / y2)
        println(name, ":")
        display(v)
    end

    mmap = NDAList([x1, x2])
    nmap = [4.2 + 0.3im, 1 / 3.0 + sqrt(2.0) * im, sin(0.7) + cos(0.4) * im]

    println("Composition of DA vectors with complex numbers.\n")
    foreach(println, compose(mmap, nmap))

    cnmap = CNDAList([y1, y2, y1 * y2])
    println("Composition of DA vectors with complex DA vectors.\n")
    foreach(display, compose(mmap, cnmap))

    cmmap = CNDAList([CNDA(x1, exp(x1)), CNDA(x2, exp(x2))])
    println("Composition of complex DA vectors with complex DA vectors.\n")
    foreach(display, compose(cmmap, cnmap))

    push!(mmap, x1 + 0.33 * x2)
    println("Composition of complex DA vectors with DA vectors.\n")
    foreach(display, compose(cmmap, mmap))
end

init!(4, 3, 400)
main()
clear!()
