# T3.2, T3.3: NDAList and the NDA algorithms, ported from python/tests/test_algorithms.py
# (Python variable 0 is Julia variable 1).

@testset "algorithms" begin
    init!(4, 3, 400)
    x = 1.0 + davar(1) + 2.0 * davar(2) + 5.0 * davar(3)
    same(a, b) = coeffs(a) == coeffs(b)

    @testset "NDAList" begin
        l = NDAList([x, exp(x)])
        push!(l, sin(x))
        @test length(l) == 3 && l isa AbstractVector{NDA}
        @test same(l[2], exp(x))
        l[1] = l[1] + 1.0                                   # getindex copy, setindex! copy
        @test con(l[1]) == con(x) + 1.0
        @test_throws BoundsError l[4]
        @test_throws BoundsError (l[0] = x)
        n = current_env().count
        l[2] = x
        @test current_env().count == n                      # set copies into the element's slot
        @test same(l[2], x)
        @test length(NDAList()) == 0
        @test same(collect(l)[3], sin(x))
    end

    @testset "compose accepts NDAList and Vector{NDA}" begin
        for as in (NDAList, identity)
            ly = compose(as([x, exp(x), sinh(x)]), as([sin(x), cos(x), tan(x)]))
            @test ly isa NDAList && length(ly) == 3
            for i in 1:3
                @test compare_da_with_file("da_composition_$(i-1).txt", ly[i], 1e-14)
            end
        end
    end

    @testset "der and integ" begin
        y = log(x)
        @test compare_da_with_file("da_der.txt", der(y, 2), 1e-14)
        @test compare_da_with_file("da_int.txt", integ(y, 2), 1e-14)
        @test_throws BoundsError der(y, 4)
        @test_throws BoundsError der(y, 0)
        @test_throws BoundsError integ(y, 4)
    end

    @testset "substitute" begin
        @test compare_da_with_file("substitute_number.txt", substitute(exp(x), 1, 1.0), 1e-14)
        @test compare_da_with_file("substitute_number.txt", substitute(exp(x), 1, 1), 1e-14)
        @test compare_da_with_file("substitute_da_vector.txt", substitute(exp(x), 1, x), 1e-14)
        z = substitute(exp(x), [1, 2], [sin(x), cos(x)])
        @test compare_da_with_file("substitute_multiple_da_vectors.txt", z, 1e-14)
        @test_throws ArgumentError substitute(exp(x), [1, 1], [sin(x), cos(x)])
        @test_throws BoundsError substitute(exp(x), [1, 4], [sin(x), cos(x)])
        @test_throws ArgumentError substitute(exp(x), [1], [sin(x), cos(x)])
        ly = substitute(NDAList([x, exp(x), sinh(x)]), [1, 2], [sin(x), cos(x)])
        @test ly isa NDAList
        for i in 1:3
            @test compare_da_with_file("bunch_substitution_$(i-1).txt", ly[i], 1e-14)
        end
    end

    @testset "compose with numbers" begin
        m = NDAList([x, exp(x), sinh(x)])
        y = compose(m, [0.1, -0.2, 0.3])
        @test y isa Vector{Float64}
        @test y[1] ≈ 1.0 + 0.1 - 0.4 + 1.5 rtol = 1e-15
        @test y[2] ≈ sum(1.2^n / factorial(n) for n in 0:4) * ℯ
        zc = compose(m, [0.1im, 0.0, 0.0])
        @test zc isa Vector{ComplexF64} && zc[1] ≈ 1.0 + 0.1im
        @test_throws ArgumentError compose(m, [0.1, 0.2])
        @test_throws ArgumentError compose(m, [x, x])
    end

    @testset "inv_map" begin
        init!(5, 2, 1000)
        a, b = davar(1), davar(2)
        m = NDAList([2.0 * a + 0.3 * b + 0.1 * a * a + 0.05 * a * b,
                     -0.4 * a + 1.5 * b + 0.2 * b * b])
        inv = inv_map(m)
        @test inv isa NDAList && length(inv) == 2
        out = compose(m, inv)
        for i in 1:2
            @test norm(out[i] - davar(i)) < 1e-12
        end
        @test_throws ArgumentError inv_map(NDAList([m[1] + 1.0, m[2]]))
        @test_throws ArgumentError inv_map(m, 3)
    end

    @testset "singular map (T3.3)" begin
        init!(5, 2, 1000)
        a, b = davar(1), davar(2)
        zero_row = NDAList([2.0 * a + a * b, 0.3 * a * a])      # linear part row 2 is zero
        @test_throws ArgumentError inv_map(zero_row)
        rank1 = [a + 2.0 * b + a * a, 2.0 * a + 4.0 * b]        # rows linearly dependent
        err = try inv_map(rank1) catch e; e end
        @test err isa ArgumentError && occursin("singular", err.msg)
    end

    init!(4, 3, 400)
    x = 1.0 + davar(1) + 2.0 * davar(2) + 5.0 * davar(3)

    @testset "evaluate_map and call" begin
        m = NDAList([x, exp(x), sinh(x)])
        pts = [sin(0.37 * (3p + j)) for j in 1:3, p in 0:999]
        out = evaluate_map(m, pts)
        @test size(out) == (3, 1000)
        loop = reduce(hcat, [compose(m, pts[:, p]) for p in 1:1000])
        @test out == loop
        @test m[2](pts[:, 7]) == out[2, 7]
        @test evaluate_map([x, exp(x)], pts)[2, 9] == out[2, 9]
        @test_throws DimensionMismatch evaluate_map(m, zeros(2, 4))
        @test_throws ArgumentError x([0.1, 0.2])
    end

    @testset "exponents" begin
        e = exponents()
        env = current_env()
        @test e isa Matrix{Int32} && size(e) == (env.nvars, env.full_length)
        v = NDA(collect(1.0:env.full_length))                 # dense
        exps = Vector{Cint}(undef, env.nvars)
        c = Ref{Cdouble}()
        for i in 1:env.full_length
            MiraDAC.check(MiraDAC.mdac_nda_index_term(v, i - 1, exps, c))
            @test e[:, i] == exps && coeff(v, e[:, i]) == i
        end
    end

    @testset "a list of a cleared env throws EnvError" begin
        init!(4, 3, 400)
        stale = NDAList([davar(1), davar(2), davar(3)])
        clear!()
        init!(4, 3, 400)
        m = NDAList([davar(1)])
        @test_throws EnvError compose(m, stale)
        @test_throws EnvError evaluate_map(stale, zeros(3, 1))
        @test_throws EnvError stale[1]
        stale = nothing
        GC.gc()
        MiraDAC.drain!()
    end
end
