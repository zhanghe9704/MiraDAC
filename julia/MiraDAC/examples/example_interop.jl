# Port of examples/example_interop.cc: NDA * SDA -> SDA, then evaluate -> NDA.
using MiraDAC

function main()
    println("=== MiraDAC Interop Example ===")
    e = current_env()
    println("order=$(e.order)  nv=$(e.nvars)  full_len=$(e.full_length)\n")

    # 1. Numerical DA vector: nx = x (the first base variable).
    nx = davar(1)
    println("NDA: nx = davar(1)")
    display(nx)

    # 2. Symbolic DA vector: sy = alpha * y.
    alpha = SymExpr("alpha")
    sy = alpha * davar(2)
    println("\nSDA: sy = alpha * davar(2)")
    display(sy)

    # 3. Mixed multiplication: NDA * SDA -> SDA (the NDA is promoted).
    result_sda = nx * sy
    println("\nSDA result = nx * sy  (NDA * SDA -> SDA, via promote):")
    display(result_sda)

    # 4. Evaluate: substitute alpha = 3.0, get an NDA.
    result_nda = evaluate(result_sda, Dict(alpha => 3.0))
    println("\nNDA result = evaluate(result_sda, Dict(alpha => 3.0)):")
    display(result_nda)

    # The coefficient of x*y should be 3.0.
    c = coeff(result_nda, [1, 1])
    println("\nCoefficient of x*y = $c  (expected 3.0)")
    if abs(c - 3.0) < 1e-12
        println("PASSED: coefficient matches 3.0")
        return 0
    end
    println("FAILED: coefficient mismatch!")
    return 1
end

init!(4, 2, 500; table=true)
status = main()
clear!()
exit(status)
