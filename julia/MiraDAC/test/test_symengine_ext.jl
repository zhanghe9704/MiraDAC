# T5.5: the SymEngine.jl extension (expressions cross as text). `import`, not `using`: SymEngine
# exports names that MiraDAC exports too (`coeff`).
import SymEngine

@testset "SymEngine.jl extension" begin
    init!(4, 3, 400)
    @test Base.get_extension(MiraDAC, :MiraDACSymEngineExt) !== nothing
    B = SymEngine.Basic
    for (jl, cpp) in (("x^2/3 + sin(x)*exp(y)", "x**2/3 + sin(x)*exp(y)"),
                      ("2^(1/3) + pi + exp(1) + 1.5*x", "2**(1/3) + pi + E + 1.5*x"),
                      ("(x + y)^(-2)", "(x + y)**(-2)"))
        @test SymExpr(B(jl)) == SymExpr(cpp)
        @test B(SymExpr(cpp)) == B(jl)
    end
    x, y = SymEngine.symbols("x y")
    @test SymExpr(x) == SymExpr("x")
    s = SDA(x)
    @test s isa SDA && con(s) == SymExpr("x")
    t = SymExpr(x) * sdavar(1) + SymExpr(y)
    @test coeffs(evaluate(t, Dict(x => 2.0, y => 0.5))) == [0.5, 2.0]
    @test coeffs(evaluate(t, Dict{B,Float64}(x => 2.0, y => 0.5))) == [0.5, 2.0]
    z = evaluate(CSDA(t, SymExpr(y) * sdavar(2)), Dict(x => 2.0, y => 0.5))   # CSDA (T6.2)
    @test z isa CNDA && coeffs(real(z)) == [0.5, 2.0] && coeffs(imag(z)) == [0.0, 0.0, 0.5]
    # SymEngine.jl's own arithmetic, next to MiraDAC's.
    @test SymEngine.expand((x + y)^2) == x^2 + 2 * x * y + y^2
    @test B(MiraDAC.expand(SymExpr((x + y)^2))) == SymEngine.expand((x + y)^2)
end
