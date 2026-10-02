# Port of examples/example_complex_symbolic.cc: symbolic complex DA (CSDA).
#
# CSDA keeps the imaginary unit outside the coefficients, so the real and imaginary parts stay
# separable and evaluate() returns a CNDA.
using MiraDAC

# C++ compare_cd_vectors: every coefficient of both parts equal up to a relative eps.
function same(a::CNDA, b::CNDA, eps)
    near(u, v) = (n = max(length(u), length(v));
                   all(isapprox.(get.(Ref(u), 1:n, 0.0), get.(Ref(v), 1:n, 0.0); rtol=eps)))
    return near(coeffs(real(a)), coeffs(real(b))) && near(coeffs(imag(a)), coeffs(imag(b)))
end

function main()
    a, b = dasymbols("a b")

    # z = (1.5 + a*x1 + 0.1*x2) + i*(0.5 + b*x1 - 0.05*x2). The constant parts are numeric so
    # that exp, sqrt, ... can branch on them.
    re = 1.5 + a * sdavar(1) + 0.1 * sdavar(2)
    ip = 0.5 + b * sdavar(1) - 0.05 * sdavar(2)
    z = CSDA(re, ip)

    println("=== Symbolic Complex DA example ===\n")
    println("z = (1.5 + a*x1 + 0.1*x2) + i*(0.5 + b*x1 - 0.05*x2)\n")

    z2 = z * z
    println("z^2 real constant part: ", Float64(con(real(z2))))  # 1.5^2 - 0.5^2 = 2.0
    println("z^2 imag constant part: ", Float64(con(imag(z2))))  # 2*1.5*0.5 = 1.5

    ez = exp(z)
    println("\nexp(z) real constant part: ", Float64(con(real(ez))),
            " (expected ~", exp(1.5) * cos(0.5), ")")
    println("exp(z) imag constant part: ", Float64(con(imag(ez))),
            " (expected ~", exp(1.5) * sin(0.5), ")")

    # Evaluate at a = 0.3, b = 0.2 and compare with the numeric computation.
    vals = Dict(a => 0.3, b => 0.2)
    n = CNDA(1.5 + 0.3 * davar(1) + 0.1 * davar(2), 0.5 + 0.2 * davar(1) - 0.05 * davar(2))
    ok = same(evaluate(ez, vals), exp(n), 1e-12)
    println("\nEvaluation cross-check (exp): ", ok ? "PASS" : "FAIL")
    ok2 = same(evaluate(sqrt(z), vals), sqrt(n), 1e-12)
    println("Evaluation cross-check (sqrt): ", ok2 ? "PASS" : "FAIL")

    (ok && ok2) || return 1
    println("\nAll CSDA operations completed successfully.")
    return 0
end

init!(3, 2, 300)
status = main()
clear!()
exit(status)
