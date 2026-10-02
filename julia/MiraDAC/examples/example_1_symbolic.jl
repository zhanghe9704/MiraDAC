# Port of examples/example_1_symbolic.cc: SDA fundamental operations.
using MiraDAC

function main()
    x, y, z = dasymbols("x y z")
    b = davar

    # An Integer stays an exact integer coefficient, a Float64 becomes a float.
    da1 = 1 + (1 + x) * b(1) + y * b(2) + (z - 0.5) * b(3)
    println("Symbolic DA vector 1:")
    display(da1)

    da2 = 3.3 + (0.5 + x) * b(1) + y * y * b(2) + (x + z + 1.1) * b(3)
    println("Symbolic DA vector 2:")
    display(da2)

    println("Summation: vec 1 + vec 2:")
    display(da1 + da2)

    println("Multiplication: vec 1 * vec 2:")
    display(da1 * da2)

    println("Inverse of vec 1:")
    display(1 / da1)

    println("Square root of vec 1:")
    display(sqrt(da1))
end

init!(4, 3, 400; table=true)
main()
clear!()
