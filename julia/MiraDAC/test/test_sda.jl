# T5.4: SDA, ported from python/tests/test_sda.py (Python variable 0 is Julia variable 1). The
# Python tests of setters (set_element, con =) have no counterpart: the Julia API has none (A.6).

scoeffs(v) = string.(coeffs(v))

@testset "SDA" begin
    init!(4, 3, 400)
    a, b = dasymbols("a b")

    @testset "constructors" begin
        @test scoeffs(SDA()) == ["0"]
        @test con(SDA(1.5)) == 1.5
        @test sprint(show, con(SDA(2))) == "SymExpr(\"2\")"      # exact integer, not 2.0
        @test con(SDA(a)) == a
        @test_throws MethodError SDA("a")
    end

    @testset "promote_sda and sdavar" begin
        x = 1.5 + 2 * davar(1) + 0.25 * davar(3)
        s = promote_sda(x)
        @test s isa SDA
        @test scoeffs(s) == ["1.5", "2", "0", "0.25"]
        @test scoeffs(sdavar(2)) == ["0", "0", "1"]
        @test_throws BoundsError sdavar(4)
        @test_throws BoundsError sdavar(0)
    end

    @testset "inspection" begin
        s = a * sdavar(1) + b * sdavar(2) * sdavar(2) + 1
        @test con(s) == 1
        @test nterms(s) == 3
        @test coeffs(s)[2] == a
        @test coeff(s, [0, 2, 0]) == b
        @test coeff(s, [0, 2]) == b
        @test_throws ArgumentError coeff(s, [-1, 0, 0])
        @test sprint(show, s) == "SDA(order=4, nvars=3, nonzero=3)"
        full = sprint(show, MIME"text/plain"(), s)
        @test full == MiraDAC.to_string(s) && occursin("a", full) && occursin("b", full)
        @test MiraDAC.env(s) == current_env()
    end

    @testset "iszero and copy" begin
        @test iszero(a * sdavar(1) - a * sdavar(1))
        t = a * sdavar(1)
        @test !iszero(t)
        for c in (copy(t), deepcopy(t))
            @test c.ptr != t.ptr && coeffs(c) == coeffs(t)
        end
    end

    @testset "binary operators" begin
        s = a + sdavar(1)                       # [a, 1]
        t = SDA(2) * sdavar(1)                  # [0, 2]
        x = davar(1)                            # NDA [0, 1]
        for other in (t, 0.5, 3, a, x)
            @test s + other isa SDA && other + s isa SDA
        end
        @test scoeffs(s + 3) == ["3 + a", "1"]
        @test scoeffs(3 + s) == ["3 + a", "1"]
        @test scoeffs(s - 1) == ["-1 + a", "1"]
        @test scoeffs(1 - s) == ["1 - a", "-1"]
        @test scoeffs(s * 2) == ["2*a", "2"]
        @test scoeffs(2 * s) == ["2*a", "2"]
        @test scoeffs(s / 2) == ["(1/2)*a", "1/2"]
        @test scoeffs(s + 0.5) == ["0.5 + a", "1"]
        @test scoeffs(s * 0.5) == ["0.5*a", "0.5"]
        @test scoeffs(s * a) == ["a**2", "a"]
        @test scoeffs(a * s) == ["a**2", "a"]
        @test scoeffs(s - a) == ["0", "1"]
        @test scoeffs(a - s) == ["0", "-1"]
        @test scoeffs(s / a) == ["1", "a**(-1)"]
        @test scoeffs(s + t) == ["a", "3"]
        @test scoeffs(s - t) == ["a", "-1"]
        @test scoeffs(s + x) == ["a", "2"]
        @test scoeffs(x + s) == ["a", "2"]
        @test scoeffs(s - x) == ["a", "0"]
        @test scoeffs(x - s) == ["-a", "0"]
        @test scoeffs(MiraDAC.expand(s * x)) == ["0", "a", "0", "0", "1"]
        @test scoeffs(-s) == ["-1.0*a", "-1.0"]
        @test scoeffs(s^2) == scoeffs(s * s)
        one = SDA(1) + sdavar(1)
        @test scoeffs(MiraDAC.expand(1 / one)) == scoeffs(MiraDAC.expand(one^-1))
        @test scoeffs(x / one) == scoeffs(promote_sda(x) / one)
        @test s^0.5 isa SDA
        @test_throws ArgumentError s / 0
        @test_throws ArgumentError s / 0.0
        @test_throws ArgumentError s / SymExpr(0)
        @test_throws MethodError s + "a"
    end

    @testset "NDA with SymExpr gives SDA" begin
        x = davar(1)
        @test scoeffs(a + x) == ["a", "1"]
        @test scoeffs(x + a) == ["a", "1"]
        @test scoeffs(x - a) == ["-a", "1"]
        @test scoeffs(a - x) == ["a", "-1"]
        @test scoeffs(x * a) == ["0", "a"]
        @test scoeffs(a * x) == ["0", "a"]
        @test scoeffs(x / a) == ["0", "a**(-1)"]
        @test a / (1 + x) isa SDA
    end

    @testset "in-place forms" begin
        base = a + sdavar(1)
        for (f!, op) in ((add!, +), (sub!, -), (mul!, *), (div!, /))
            for other in (SDA(2) + sdavar(2), 2.0, 2, a, 2 + davar(2))
                want = scoeffs(MiraDAC.expand(op(base, other)))
                s = copy(base)
                p = s.ptr
                @test f!(s, s, other) === s && s.ptr == p
                @test scoeffs(MiraDAC.expand(s)) == want
                o = SDA()
                f!(o, base, other)
                @test scoeffs(MiraDAC.expand(o)) == want
                f!(o, other, base)
                @test scoeffs(MiraDAC.expand(o)) == scoeffs(MiraDAC.expand(op(other, base)))
            end
        end
        o = SDA()
        add!(o, davar(1), a)                    # NDA with SymExpr
        @test scoeffs(o) == ["a", "1"]
        exp!(o, base)
        @test scoeffs(o) == scoeffs(exp(base))
        @test_throws ArgumentError div!(o, base, 0.0)
    end

    @testset "env checks" begin
        s = sdavar(1)
        clear!()
        init!(4, 3, 400)
        @test_throws EnvError s + 1
        @test_throws EnvError sdavar(1) + s
        @test sprint(show, s) == "SDA(cleared env)"
    end

    @testset "functions" begin
        s = a + sdavar(1)
        for f in (sqrt, exp, log, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, erf)
            @test f(s) isa SDA
        end
        @test scoeffs(exp(s))[1:2] == ["exp(a)", "1.0*exp(a)"]
        for f in (asinh, acosh, atanh, abs)
            @test_throws MethodError f(s)
        end
        @test exp(davar(1)) isa NDA             # the NDA methods are unchanged
    end

    @testset "der and integ" begin
        s = a * sdavar(1) * sdavar(1)
        @test scoeffs(der(s, 1)) == ["0", "2.0*a"]
        @test coeff(integ(sdavar(2) * a, 2), [0, 2, 0]) == a / 2
        @test_throws BoundsError der(s, 4)
        @test_throws BoundsError integ(s, 0)
    end

    @testset "substitute, compose and SDAList" begin
        s = a * sdavar(1) + sdavar(2)
        @test scoeffs(substitute(s, 1, 2.0)) == scoeffs(2.0 * a + sdavar(2))
        @test scoeffs(substitute(s, 2, sdavar(1))) == scoeffs(a * sdavar(1) + sdavar(1))
        @test scoeffs(substitute(s, [1, 2], SDAList([SDA(1), SDA(a)]))) == ["2*a"]
        @test scoeffs(substitute(s, [1, 2], [SDA(1), SDA(a)])) == ["2*a"]
        @test_throws ArgumentError substitute(s, [1, 1], [SDA(1), SDA(1)])
        m = [s, s * s]
        pt = [SDA(1), SDA(2), SDA(0)]
        outs = compose(m, pt)
        @test outs isa SDAList
        @test [scoeffs(v) for v in outs] == [["2 + a"], ["4 + 4*a + a**2"]]
        @test_throws ArgumentError compose(m, pt[1:2])
        sub = substitute(m, [1, 2], [SDA(1), SDA(a)])
        @test sub isa SDAList && scoeffs(sub[1]) == ["2*a"]
        vals = compose(m, [1.0, 2.0, 0.0])
        @test vals isa Vector{SymExpr}
        @test [Float64(MiraDAC.subs(e, Dict(a => 0.5))) for e in vals] == [2.5, 6.25]

        l = SDAList(m)
        @test length(l) == 2 && scoeffs(l[2]) == scoeffs(s * s)
        l[2] = s
        @test scoeffs(l[2]) == scoeffs(s)
        push!(l, sdavar(3))
        @test length(l) == 3 && MiraDAC.env(l) == current_env()
        @test_throws EnvError MiraDAC.env(SDAList())
    end

    @testset "evaluate" begin
        s = exp(a * sdavar(1) + 1)
        got = evaluate(s, Dict(a => 0.3))
        want = exp(0.3 * davar(1) + 1)
        @test got isa NDA
        @test length(coeffs(got)) == length(coeffs(want))
        @test all(abs(p - q) <= 1e-13 * max(1.0, abs(q)) for (p, q) in zip(coeffs(got), coeffs(want)))
        x = 1.5 + 2.0 * davar(1) + 0.5 * davar(2) * davar(3)
        @test coeffs(evaluate(promote_sda(x), Dict{SymExpr,Float64}())) == coeffs(x)
        @test coeffs(evaluate(a * sdavar(1) + b, Dict(a => 2, b => 3.5))) == [3.5, 2.0]
        @test_throws ArgumentError evaluate(s, Dict{SymExpr,Float64}())    # a has no value
    end

    @testset "simplify, expand, subs" begin
        s = (a + 1)^2 * sdavar(1) + con(sin(SDA(a))^2) + con(cos(SDA(a))^2)
        e = MiraDAC.expand(s)
        @test e.ptr != s.ptr
        @test coeffs(e)[2] == a^2 + 2 * a + 1
        @test coeffs(s)[2] == (a + 1)^2                     # unchanged
        t = (a * 2 / 2 + b - b) * sdavar(2)
        @test coeff(simplify(t), [0, 1, 0]) == a
        @test coeffs(MiraDAC.subs(s, Dict(a => 2)))[2] == 9
        @test coeffs(MiraDAC.subs(s, Dict(a => b)))[2] == (b + 1)^2
    end

    @testset "pool retry" begin
        init!(4, 3, 32)
        r = MiraDAC.POOL_RETRIES[]
        x = SDA(a)
        for _ in 1:500
            x = x * 0.5 + sdavar(1)
        end
        @test Float64(coeffs(x)[2]) ≈ 2.0
        @test MiraDAC.POOL_RETRIES[] > r
        # compose at a point takes nvars + length(m) temporary slots
        init!(3, 2, 40)
        m = SDAList([sdavar(1), sdavar(2)])
        GC.gc(); MiraDAC.drain!()
        r, k = MiraDAC.POOL_RETRIES[], 0
        while MiraDAC.POOL_RETRIES[] == r       # count the garbage that fills the pool
            SDA(1.0); k += 1
        end
        GC.gc(); MiraDAC.drain!()
        for _ in 1:k-2
            SDA(1.0)
        end
        r = MiraDAC.POOL_RETRIES[]
        @test compose(m, [1.0, 2.0]) == SymExpr[SymExpr(1.0), SymExpr(2.0)]
        @test MiraDAC.POOL_RETRIES[] > r
        init!(4, 3, 400)
    end
end
