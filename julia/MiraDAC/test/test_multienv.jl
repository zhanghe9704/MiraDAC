# T7.3: the multi-env API, ported from python/tests/test_multienv.py (Python variable 0 is
# Julia variable 1).

@testset "multienv" begin
    env = MiraDAC.env
    # A: the default env (order 2, 2 vars), current. B: a second env (order 4, 3 vars).
    function envs(f)
        init!(2, 2, 200)
        b = DAEnv(4, 3, 4000)
        try
            f(default_env(), b)
        finally
            close(b)
            clear!()
        end
    end

    @testset "creation and properties" begin
        envs() do a, b
            @test current_env() == a                     # creating b kept a current
            @test (b.order, b.max_order, b.nvars, b.full_length) == (4, 4, 3, 35)
            @test (b.poolsize, b.count, b.remain) == (4000, 0, 4000)
            @test !b.retired && a != b && length(Set([a, b, default_env()])) == 2
            @test hash(a) == hash(default_env()) && hash(a) != hash(b)
            @test sprint(show, b) == "DAEnv(order=4, nvars=3, poolsize=4000)"
            x = davar(1; env=b)
            @test b.count == 1 && env(x) == b && env(davar(1)) == a
            @test_throws ArgumentError DAEnv(2, 2, 0)
        end
    end

    @testset "current and default without env" begin
        clear!()
        @test_throws EnvError current_env()
        e = DAEnv(2, 2, 10)
        @test_throws EnvError current_env()          # creating an env does not select it
        with_env(e) do
            @test current_env() == e
        end
        @test_throws EnvError current_env()
        close(e)
        @test_throws EnvError with_env(() -> nothing, e)
        init!(2, 2, 10)
        clear!()
    end

    @testset "with_env nests and restores" begin
        envs() do a, b
            with_env(b) do
                @test current_env() == b
                with_env(a) do
                    with_env(b) do
                        @test env(davar(1)) == b
                    end
                    @test current_env() == a
                end
                with_env(b) do                       # re-entrant
                    @test current_env() == b
                end
                @test current_env() == b
            end
            @test current_env() == a
            @test with_env(() -> :r, b) === :r
            @test_throws ErrorException with_env(b) do
                error("boom")
            end
            @test current_env() == a
        end
    end

    @testset "env keyword on constructors" begin
        envs() do a, b
            made = Any[NDA(1.5; env=b), NDA(ones(4); env=b), davar(3; env=b), CNDA(1 + 2im; env=b)]
            if MiraDAC.HAS_SYMBOLIC
                x, = dasymbols("x")
                append!(made, [SDA(; env=b), SDA(2.5; env=b), SDA(2; env=b), SDA(x; env=b),
                               sdavar(3; env=b), CSDA(; env=b), CSDA(1 + 2im; env=b)])
            end
            @test all(v -> env(v) == b, made)
            @test env(NDA(1.5)) == a
            @test size(exponents(b)) == (3, 35) && size(exponents()) == (2, 6)
            @test_throws BoundsError davar(3)        # a has 2 variables, b has 3
            @test current_env() == a
        end
    end

    @testset "operations on env-B vectors while A is current" begin
        envs() do a, b
            xa, xb = davar(1), davar(1; env=b)
            @test nterms(xa * xa * xa) == 0                          # order 2
            @test coeff(xb * xb * xb, [3, 0, 0]) == 1.0              # order 4
            @test nterms(xb^5) == 0 && nterms(xb^4) == 1
            @test nterms(exp(xa)) == 3 && nterms(exp(xb)) == 5

            x = 0.3 + 0.5 * davar(1; env=b) + davar(2; env=b) - 0.2 * davar(3; env=b)
            m = [2.0 * davar(1; env=b) + 0.1 * davar(2; env=b)^2, davar(2; env=b),
                 davar(3; env=b) + 0.1 * davar(1; env=b) * davar(3; env=b)]
            out = NDA(0.0; env=b)
            cases = [
                () -> coeffs(x * x + x / 2.0 - 1.0),
                () -> coeffs(exp(x)), () -> coeffs(sin(x)), () -> coeffs(MiraDAC.erf(x)),
                () -> coeffs(mul!(out, x, x)), () -> coeffs(exp!(out, x)),
                () -> coeffs(der(sin(x), 2)), () -> coeffs(integ(sin(x), 2)),
                () -> coeffs(substitute(sin(x), 1, 0.5)),
                () -> coeffs.(compose([x, sin(x), cos(x)], [cos(x), x, sin(x)])),
                () -> coeffs.(inv_map(m)),
                () -> evaluate_map([x, sin(x)], fill(0.1, 3, 4)),
                () -> (c = CNDA(x, 2.0 * x) * CNDA(x); (coeffs(real(c)), coeffs(imag(c)))),
                () -> (c = exp(CNDA(x, x)); (coeffs(real(c)), coeffs(imag(c)))),
            ]
            if MiraDAC.HAS_SYMBOLIC
                p, q = dasymbols("p q")
                s = p + 0.5 * sdavar(1; env=b) + q * sdavar(2; env=b) * sdavar(3; env=b)
                append!(cases, [
                    () -> string.(coeffs(s * s + 2.0)),
                    () -> string.(coeffs(exp(s))),
                    () -> string.(coeffs(der(s, 1))),
                    () -> coeffs(evaluate(exp(s), Dict(p => 0.3, q => -0.2))),
                    () -> string(CSDA(s, s) * CSDA(s)),
                ])
            end
            for f in cases
                want = with_env(f, b)
                @test current_env() == a
                @test f() == want
            end
            @test env(exp(x)) == b && env(CNDA(x)) == b
            @test current_env() == a
        end
    end

    @testset "mixing envs throws EnvError" begin
        envs() do a, b
            xa, xb = davar(1), davar(1; env=b)
            for op in (() -> xa + xb, () -> xb * xa, () -> add!(xb, xb, xa), () -> CNDA(xa, xb),
                       () -> CNDA(xb) * CNDA(xa),
                       () -> compose([xb, xb, xb], [xa, xa, xa]))
                @test_throws EnvError op()
            end
            if MiraDAC.HAS_SYMBOLIC
                sa, sb = sdavar(1), sdavar(1; env=b)
                for op in (() -> sa + sb, () -> sb * xa, () -> CSDA(sa, sb))
                    @test_throws EnvError op()
                end
            end
        end
    end

    @testset "close with live vectors" begin
        envs() do a, b
            x = exp(davar(1; env=b))
            c = CNDA(x, x)
            close(b)
            @test b.retired && sprint(show, b) == "DAEnv(retired)"
            close(b)                                 # idempotent
            for use in (() -> x + 1.0, () -> exp(x), () -> con(x), () -> c * c,
                        () -> NDA(1.0; env=b), () -> b.count, () -> import_vec(a, x),
                        () -> import_vec(b, davar(1)), () -> with_env(() -> nothing, b),
                        () -> with_order(() -> nothing, b, 1))
                @test_throws EnvError use()
            end
            @test env(x) == b && current_env() == a
            if MiraDAC.HAS_SYMBOLIC
                @test_throws EnvError promote_sda(b, davar(1))
            end
            x = c = nothing
            GC.gc()                                  # finalizers of the closed env's vectors
            MiraDAC.drain!()
            @test current_env() == a
        end
    end

    @testset "close the default env" begin
        envs() do a, b
            x = davar(1)
            close(a)
            @test a.retired
            @test_throws EnvError default_env()
            @test_throws EnvError x * x
            init!(2, 2, 10)                          # envs() ends with clear!()
        end
    end

    @testset "import_vec and promote_sda" begin
        envs() do a, b
            b2 = DAEnv(4, 3, 100)
            try
                x = 1.5 + davar(1; env=b) * davar(2; env=b) - 2.0 * davar(3; env=b)
                y = import_vec(b2, x)
                @test env(y) == b2 && b2.count == 1 && current_env() == a
                z = import_vec(b, y)
                @test env(z) == b && coeffs(z) == coeffs(x)
                c = CNDA(x, 2.0 * x)
                cc = import_vec(b, import_vec(b2, c))
                @test env(cc) == b && coeffs(real(cc)) == coeffs(x)
                @test coeffs(imag(cc)) == coeffs(2.0 * x)
                @test_throws ArgumentError import_vec(a, x)          # different layout
                if MiraDAC.HAS_SYMBOLIC
                    s = promote_sda(b2, x)
                    @test env(s) == b2
                    @test string.(coeffs(s)) == string.(coeffs(promote_sda(x)))
                    t = import_vec(b, s)
                    @test env(t) == b && string.(coeffs(t)) == string.(coeffs(s))
                    cs = promote_sda(b2, c)
                    @test cs isa CSDA && env(cs) == b2 && env(real(cs)) == b2
                    cs2 = import_vec(b, cs)
                    @test env(cs2) == b
                    @test string.(coeffs(imag(cs2))) == string.(coeffs(imag(cs)))
                    @test_throws ArgumentError promote_sda(a, x)
                    @test_throws ArgumentError import_vec(a, s)
                end
            finally
                close(b2)
            end
        end
    end

    @testset "with_order affects one env only" begin
        envs() do a, b
            xb = davar(1; env=b)
            with_order(b, 2) do
                @test (b.order, a.order, current_env().order) == (2, 2, 2)
                @test nterms(xb * xb * xb) == 0
                with_order(b, 1) do
                    @test b.order == 1 && nterms(xb * xb) == 0
                end
                @test b.order == 2
                with_order(1) do                     # the current env, a
                    @test (a.order, b.order) == (1, 2)
                end
                @test a.order == 2
            end
            @test b.order == 4 && coeff(xb * xb * xb, [3, 0, 0]) == 1.0
            @test_throws ArgumentError with_order(() -> nothing, b, 5)
            @test_throws KeyError with_order(b, 3) do
                throw(KeyError(:k))
            end
            @test b.order == 4
        end
    end

    @testset "inv_map in a non-default env" begin
        function inv_coeffs()
            x, y = davar(1), davar(2)
            m = [2.0 * x + 0.3 * y + 0.1 * x * x + 0.05 * x * y, -0.4 * x + 1.5 * y + 0.2 * y * y]
            return coeffs.(inv_map(m))
        end
        init!(5, 2, 1000)
        want = inv_coeffs()
        e = DAEnv(5, 2, 1000)
        got = with_env(inv_coeffs, e)
        clear!()
        again = with_env(inv_coeffs, e)               # also with no default env
        close(e)
        @test got == want && again == want
    end

    @testset "default env without user vectors is retired, not freed" begin
        init!(2, 2, 10)
        a = default_env()
        clear!()
        @test a.retired
        init!(3, 3, 20)
        @test default_env() != a && a.retired
        b = default_env()
        close(b)
        @test b.retired
        clear!()
    end

    # A 600 GB pool fails at once only where the kernel refuses to overcommit that much
    # (Linux); macOS hands out the address space lazily, and filling it gets the process killed.
    Sys.islinux() && @testset "failed env construction keeps the current env" begin
        envs() do a, b
            @test_throws ErrorException DAEnv(4, 3, 2^31)    # std::bad_alloc, MDAC_ERR_RUNTIME
            @test current_env() == a
        end
    end

    if MiraDAC.HAS_SYMBOLIC
        @testset "erf of an SDA in a non-default table env" begin
            # The last case of test/test_symbolic.cc.
            sx, sy = dasymbols("x y")
            function erf_coeffs()
                s = MiraDAC.erf(sx + sdavar(1) * sdavar(1) + sy * sdavar(2))
                return coeffs(evaluate(s, Dict(sx => 0.3, sy => -0.7)))
            end
            init!(5, 2, 400; table=true)
            want = erf_coeffs()
            e = DAEnv(5, 2, 400; table=true)
            got = with_env(erf_coeffs, e)
            close(e)
            clear!()
            @test length(got) == length(want)
            @test all(abs.(got .- want) .< 1e-13)
            @test want[1] ≈ 0.3286267594591274                  # erf(0.3)
        end
    end

    init!(4, 3, 400)
end
