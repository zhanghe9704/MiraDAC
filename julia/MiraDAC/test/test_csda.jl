# T6.2: CSDA, ported from python/tests/test_csda.py (the port of test_symbolic_cd.cc and the
# checks of examples/example_complex_symbolic.cc; Python variable 0 is Julia variable 1). The
# Python tests of part setters (`z.imag = ...`) have no counterpart: the Julia API has none (A.6);
# nor has the list form of `evaluate`, nor CSDA operators with a `SymEngine.Basic`.

make_csda(c0r, b1r, b2r, c0i, b1i, b2i) =
    CSDA(SDA(c0r) + b1r * sdavar(1) + b2r * sdavar(2), SDA(c0i) + b1i * sdavar(1) + b2i * sdavar(2))
make_cnda(c0r, b1r, b2r, c0i, b1i, b2i) =
    CNDA(NDA(c0r) + b1r * davar(1) + b2r * davar(2), NDA(c0i) + b1i * davar(1) + b2i * davar(2))

@testset "CSDA" begin
    init!(3, 2, 800)

    # ---- port of test_symbolic_cd.cc ----

    @testset "structural check" begin
        (b,) = dasymbols("b")
        z = make_csda(1.0, b, 0.0, 0.5, 0.2, 0.0)
        @test Float64(con(real(z))) ≈ 1.0
        @test Float64(con(imag(z))) ≈ 0.5
        d = z - z
        @test Float64(con(real(d))) == 0.0 && Float64(con(imag(d))) == 0.0
    end

    @testset "arithmetic cross-check" begin
        syms = dasymbols("b1 b2 c1 c2 e1 e2 f1 f2")
        b1, b2, c1, c2, e1, e2, f1, f2 = syms
        z1 = make_csda(0.8, b1, b2, 0.3, c1, c2)
        z2 = make_csda(0.5, e1, e2, 0.6, f1, f2)
        v = [0.3, -0.1, 0.2, 0.4, -0.15, 0.05, 0.1, -0.2]
        vals = Dict(zip(syms, v))
        n1 = make_cnda(0.8, v[1], v[2], 0.3, v[3], v[4])
        n2 = make_cnda(0.5, v[5], v[6], 0.6, v[7], v[8])
        c = 0.5 - 0.3im
        for (rs, rn) in ((z1 + z2, n1 + n2), (z1 - z2, n1 - n2), (z1 * z2, n1 * n2),
                         (z1 / z2, n1 / n2), (z1 + c, n1 + c), (z1 * 2.5, n1 * 2.5),
                         (1.0 / z2, 1.0 / n2), ((-1.0) * z1, (-1.0) * n1))
            @test near(evaluate(rs, vals), rn, 1e-10)
        end
    end

    @testset "function cross-check: $fn" for (fn, cons) in (
            (exp, (0.5, 0.3)), (sqrt, (0.5, 0.3)), (log, (0.5, 0.3)), (asin, (0.3, 0.2)),
            (acos, (0.3, 0.2)), (atan, (0.3, 0.2)), (asinh, (0.5, 0.3)), (atanh, (0.5, 0.3)),
            (acosh, (1.5, 0.2)))
        b1, c1 = dasymbols("b1 c1")
        z = make_csda(cons[1], b1, 0.0, cons[2], c1, 0.0)
        n = make_cnda(cons[1], 0.08, 0.0, cons[2], 0.05, 0.0)
        @test near(evaluate(fn(z), Dict(b1 => 0.08, c1 => 0.05)), fn(n), 1e-9)
    end

    @testset "exp, sqrt, log with two symbols" begin
        b1, b2, c1, c2 = dasymbols("b1 b2 c1 c2")
        z = make_csda(0.5, b1, b2, 0.3, c1, c2)
        vals = Dict(b1 => 0.10, b2 => 0.05, c1 => 0.08, c2 => -0.04)
        n = make_cnda(0.5, 0.10, 0.05, 0.3, 0.08, -0.04)
        for f in (exp, sqrt, log)
            @test near(evaluate(f(z), vals), f(n), 1e-9)
        end
    end

    @testset "pow cross-check" begin
        b1, c1 = dasymbols("b1 c1")
        z = make_csda(0.5, b1, 0.0, 0.3, c1, 0.0)
        vals = Dict(b1 => 0.08, c1 => 0.05)
        n = make_cnda(0.5, 0.08, 0.0, 0.3, 0.05, 0.0)
        @test near(evaluate(z^3, vals), n^3, 1e-9)
        @test near(evaluate(z^2.5, vals), n^2.5, 1e-9)
        @test near(evaluate(z^-1, vals), n^-1, 1e-9)
    end

    @testset "abs is the symbolic modulus" begin
        b1, c1 = dasymbols("b1 c1")
        z = make_csda(0.6, b1, 0.0, 0.4, c1, 0.0)
        n = make_cnda(0.6, 0.1, 0.0, 0.4, 0.08, 0.0)
        r = abs(z)
        @test r isa SDA
        @test near(evaluate(r, Dict(b1 => 0.1, c1 => 0.08)), sqrt(real(n) * real(n) + imag(n) * imag(n)), 1e-10)
    end

    @testset "promote_sda and evaluate round trip" begin
        orig = CNDA(1.5 + 0.3 * davar(1) - 0.1 * davar(2), 0.7 + 0.2 * davar(1) - 0.05 * davar(2))
        p = promote_sda(orig)
        @test p isa CSDA
        @test near(evaluate(p, Dict{SymExpr,Float64}()), orig, 1e-14)
    end

    @testset "cd_composition cross-check" begin
        b1, c1 = dasymbols("b1 c1")
        z0 = make_csda(1.2, b1, 0.0, 0.4, c1, 0.0)
        z1 = make_csda(0.5, 0.0, 0.0, 0.2, 0.0, 0.0)
        v0re, v0im = 0.7 + 0.1 * davar(1), 0.05 * davar(1)
        v1re, v1im = 0.3 + 0.05 * davar(2), 0.02 * davar(2)
        cvec = CSDAList([CSDA(promote_sda(v0re), promote_sda(v0im)),
                         CSDA(promote_sda(v1re), promote_sda(v1im))])
        out = compose([z0, z1], cvec)
        @test out isa CSDAList && length(out) == 2
        n0 = CNDA(1.2 + 0.3 * davar(1), 0.4 + 0.15 * davar(1))
        n1 = CNDA(NDA(0.5), NDA(0.2))
        nout = compose([n0, n1], [CNDA(v0re, v0im), CNDA(v1re, v1im)])
        vals = Dict(b1 => 0.3, c1 => 0.15)
        for (s, n) in zip(out, nout)
            @test near(evaluate(s, vals), n, 1e-9)
        end
        # a CSDA map at SDA arguments
        sv = [promote_sda(v0re), promote_sda(v1re)]
        sout = compose([z0, z1], sv)
        nsout = compose([n0, n1], [v0re, v1re])
        for (s, n) in zip(sout, nsout)
            @test near(evaluate(s, vals), n, 1e-9)
        end
    end

    @testset "cd_composition of an SDA map at CSDA arguments" begin
        b1, c1 = dasymbols("b1 c1")
        ivec = SDAList([1.1 + b1 * sdavar(1), 0.5 + c1 * sdavar(2)])
        cv0re, cv1re = 0.5 + 0.1 * davar(1), NDA(0.3)
        cv0im, cv1im = NDA(0.2), 0.05 * davar(2)
        cmap = [CSDA(promote_sda(cv0re), promote_sda(cv0im)), CSDA(promote_sda(cv1re), promote_sda(cv1im))]
        out = compose(ivec, cmap)
        nout = compose([1.1 + 0.3 * davar(1), 0.5 + 0.2 * davar(2)], [CNDA(cv0re, cv0im), CNDA(cv1re, cv1im)])
        vals = Dict(b1 => 0.3, c1 => 0.2)
        for (s, n) in zip(out, nout)
            @test near(evaluate(s, vals), n, 1e-9)
        end
    end

    # ---- checks of example_complex_symbolic.cc ----

    @testset "example_complex_symbolic" begin
        a, b = dasymbols("a b")
        z = CSDA(SDA(1.5) + a * sdavar(1) + SDA(0.1) * sdavar(2),
                 SDA(0.5) + b * sdavar(1) + SDA(-0.05) * sdavar(2))
        z2 = z * z
        @test Float64(con(real(z2))) ≈ 2.0
        @test Float64(con(imag(z2))) ≈ 1.5
        ez = exp(z)
        @test Float64(con(real(ez))) ≈ exp(1.5) * cos(0.5) rtol = 1e-14
        @test Float64(con(imag(ez))) ≈ exp(1.5) * sin(0.5) rtol = 1e-14
        vals = Dict(a => 0.3, b => 0.2)
        n = CNDA(1.5 + 0.3 * davar(1) + 0.1 * davar(2), 0.5 + 0.2 * davar(1) - 0.05 * davar(2))
        @test near(evaluate(ez, vals), exp(n), 1e-12)
        @test near(evaluate(sqrt(z), vals), sqrt(n), 1e-12)
    end

    # ---- the type ----

    @testset "constructors, parts and show" begin
        (a,) = dasymbols("a")
        re = a * sdavar(1) + 1
        z = CSDA(re)
        @test scoeffs(real(z)) == scoeffs(re) && scoeffs(imag(z)) == ["0"]
        @test scoeffs(real(CSDA())) == ["0"] && scoeffs(imag(CSDA())) == ["0"]
        w = CSDA(1.5 - 2im)
        @test con(real(w)) == 1.5 && con(imag(w)) == -2.0
        z = CSDA(re, 2 * re)
        @test scoeffs(imag(z)) == scoeffs(2 * re)
        r = real(z)                                       # a copy
        add!(r, r, 1)
        @test con(real(z)) == 1
        @test sprint(show, z) == "CSDA(order=3, nvars=2, nonzero=(2, 2))"
        full = sprint(show, MIME"text/plain"(), z)
        @test full == MiraDAC.to_string(z)
        @test occursin("Real part", full) && occursin("Imaginary part", full) && occursin("a", full)
        for c in (copy(z), deepcopy(z))
            @test c.ptr != z.ptr && scoeffs(real(c)) == scoeffs(real(z)) && scoeffs(imag(c)) == scoeffs(imag(z))
        end
        @test MiraDAC.env(z) == current_env()
        @test_throws MethodError CSDA(davar(1))
    end

    @testset "scalar operators match CNDA" begin
        x = 0.4 + 0.3 * davar(1) - 0.2 * davar(2)
        y = 0.1 - 0.5 * davar(1) * davar(2)
        n = CNDA(x, y)
        z = promote_sda(n)
        w = 0.5 - 2.0im
        for f in (u -> u + w, u -> w - u, u -> u * w, u -> w / u, u -> u / w, u -> 2 + u,
                  u -> 2.0 - u, u -> u * 3, u -> 2.0 / u, u -> -u, u -> u^0.5)
            r = f(z)
            @test r isa CSDA
            @test near(evaluate(r, Dict{SymExpr,Float64}()), f(n), 1e-13)
        end
        @test_throws ArgumentError z / 0.0
    end

    @testset "SDA with a complex gives a CSDA" begin
        x = 0.4 + 0.3 * davar(1)
        s = promote_sda(x)
        w = 0.5 - 2.0im
        for (r, n) in ((s + w, x + w), (w + s, w + x), (s - w, x - w), (w - s, w - x),
                       (s * w, x * w), (w * s, w * x), (s / w, x / w), (w / s, w / x))
            @test r isa CSDA
            @test near(evaluate(r, Dict{SymExpr,Float64}()), n, 1e-13)
        end
        @test s + 2 isa SDA && s * 2.0 isa SDA
    end

    @testset "SDA and SymExpr operands" begin
        b1, b2, c1, c2, e1, k = dasymbols("b1 b2 c1 c2 e1 k")
        z = make_csda(0.8, b1, b2, 0.3, c1, c2)
        s = SDA(0.5) + e1 * sdavar(1)
        vals = Dict(b1 => 0.3, b2 => -0.1, c1 => 0.2, c2 => 0.4, e1 => -0.15, k => 1.7)
        zn = make_cnda(0.8, 0.3, -0.1, 0.3, 0.2, 0.4)
        sn = NDA(0.5) - 0.15 * davar(1)
        for op in (+, -, *, /), (sym, num) in ((s, sn), (k, 1.7))
            got = op(z, sym)
            @test got isa CSDA
            @test near(evaluate(got, vals), op(zn, num), 1e-12)
            @test near(evaluate(op(sym, z), vals), op(num, zn), 1e-12)
        end
    end

    @testset "in-place forms: $rhs" for rhs in ("csda", 0.5 - 2.0im, 2.0, 3, "sda", "expr")
        (a,) = dasymbols("a")
        z1 = make_csda(0.8, a, 0.1, 0.3, 0.2, a)
        z2 = make_csda(0.5, 0.3, a, 0.6, 0.1, 0.2)
        other = rhs == "csda" ? z2 : rhs == "sda" ? real(z2) : rhs == "expr" ? a + 1 : rhs
        for (op, op!) in ((+, add!), (-, sub!), (*, mul!), (/, div!))
            out = copy(z1)
            p = out.ptr
            @test op!(out, out, other) === out && out.ptr == p
            @test near(evaluate(out, Dict(a => 0.7)), evaluate(op(z1, other), Dict(a => 0.7)), 1e-13)
            @test near(evaluate(op!(out, other, z1), Dict(a => 0.7)), evaluate(op(other, z1), Dict(a => 0.7)), 1e-13)
        end
        out = CSDA()
        @test near(evaluate(exp!(out, z2), Dict(a => 0.7)), evaluate(exp(z2), Dict(a => 0.7)), 1e-13)
    end

    @testset "unsupported operands and functions" begin
        z = CSDA(sdavar(1))
        for other in (davar(1), CNDA(davar(1)))
            @test_throws MethodError z + other
            @test_throws MethodError other * z
        end
        z = CSDA(SDA(0.5) + sdavar(1))
        for f in (sin, cos, tan, sinh, cosh, tanh, MiraDAC.erf)
            @test_throws MethodError f(z)
        end
        @test sin(sdavar(1)) isa SDA
        @test exp(CNDA(davar(1))) isa CNDA
    end

    @testset "evaluate errors" begin
        a, b = dasymbols("a b")
        z = make_csda(0.5, a, 0.0, 0.3, 0.0, b)
        n = make_cnda(0.5, 0.1, 0.0, 0.3, 0.0, 0.2)
        r = evaluate(z, Dict(a => 0.1, b => 0.2))
        @test r isa CNDA && near(r, n, 1e-15)
        @test_throws ArgumentError evaluate(z, Dict(a => 0.1))
    end

    @testset "cd_composition checks" begin
        z = [CSDA(sdavar(1)), CSDA(sdavar(2))]
        v = CSDAList([CSDA(sdavar(2)), CSDA(sdavar(1))])
        @test_throws ArgumentError compose(z, v[1:1])
        @test isempty(compose(CSDA[], v))
        out = compose(z, v)                               # swaps the variables
        @test near(evaluate(out[1], Dict{SymExpr,Float64}()), CNDA(davar(2)), 0)
        @test near(evaluate(out[2], Dict{SymExpr,Float64}()), CNDA(davar(1)), 0)
    end

    @testset "CSDAList" begin
        l = CSDAList([CSDA(sdavar(1))])
        push!(l, CSDA())
        @test length(l) == 2 && l[1] isa CSDA
        l[2] = CSDA(sdavar(2))
        @test scoeffs(real(l[2])) == scoeffs(sdavar(2))
        @test_throws BoundsError l[3]
        @test MiraDAC.env(l) == current_env()
        @test_throws EnvError MiraDAC.env(CSDAList())
        @test collect(CSDAList(CSDA[])) == CSDA[]
    end

    @testset "parts must share the env" begin
        s = sdavar(1)
        z = CSDA(s)
        clear!()
        init!(3, 2, 100)
        t = sdavar(2)
        @test_throws EnvError CSDA(t, s)
        @test_throws EnvError z + 1.0
        @test sprint(show, z) == "CSDA(cleared env)"
    end

    @testset "pool exhaustion retries" begin
        init!(3, 2, 40)
        r = MiraDAC.POOL_RETRIES[]
        x = CSDA(sdavar(1))
        for _ in 1:300
            x = x * 0.5 + (1.0 + 1.0im)
        end
        @test Float64(con(real(x))) ≈ 2.0 && Float64(con(imag(x))) ≈ 2.0
        @test MiraDAC.POOL_RETRIES[] > r
    end

    init!(4, 3, 400)
end
