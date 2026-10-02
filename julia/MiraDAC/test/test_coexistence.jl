# T5.6: MiraDAC's libsymengine and SymEngine.jl's (another version) in one process, loaded in
# either order; each order runs in its own Julia process.
const COEXIST_SCRIPT = raw"""
if ARGS[1] == "SymEngine"
    import SymEngine
    using MiraDAC
else
    using MiraDAC
    import SymEngine
end
init!(4, 2, 400)
a = SymExpr("a")
s = exp(a + sdavar(1)) * (1 + davar(2))
ok = [string.(coeffs(s))[1:3] == ["exp(a)", "1.0*exp(a)", "exp(a)"],
      coeffs(evaluate(s, Dict(a => 0.5)))[1:3] ≈ [exp(0.5), exp(0.5), exp(0.5)]]
x, y = SymEngine.symbols("x y")
e = SymEngine.expand((x + y)^3)
push!(ok, e == x^3 + 3 * x^2 * y + 3 * x * y^2 + y^3,
      SymEngine.subs(e, x => 1, y => 2) == 27,
      SymExpr(e) == MiraDAC.expand(SymExpr("(x + y)**3")),
      SymEngine.Basic(con(s)) == exp(SymEngine.Basic("a")))
println(all(ok) ? "COEXIST OK" : "COEXIST FAIL $ok")
"""

@testset "SymEngine coexistence" begin
    for first in ("SymEngine", "MiraDAC")
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $COEXIST_SCRIPT $first`
        out = read(pipeline(ignorestatus(cmd); stderr=stdout), String)
        @test occursin("COEXIST OK", out) || (println(out); false)
    end
end
