# T4.2: CNDA, ported from python/tests/test_cnda.py (Python variable 0 is Julia variable 1).

# Port of da::compare_cd_with_file: the real and the imaginary table compared as by
# compare_da_with_file (test_nda.jl).
function read_cd_file(path)
    re, im = Dict{Vector{Int},Float64}(), Dict{Vector{Int},Float64}()
    reading = false
    for line in eachline(path)
        line = strip(line)
        isempty(line) && continue
        if reading
            w = split(line)
            k = parse.(Int, w[4:end-1])
            re[k], im[k] = parse(Float64, w[2]), parse(Float64, w[3])
        elseif isempty(replace(line, '-' => ""))
            reading = true
        end
    end
    return re, im
end

function compare_cd_with_file(name, c::CNDA, eps)
    re, im = read_cd_file(ref_path(name))
    return compare_terms(re, nda_terms(real(c)), eps) && compare_terms(im, nda_terms(imag(c)), eps)
end

dense(v::NDA) = (c = coeffs(v); [c; zeros(MiraDAC.env(v).full_length - length(c))])
near(a::NDA, b::NDA, tol=1e-13) = tol == 0 ? dense(a) == dense(b) :
                                  isapprox(dense(a), dense(b); rtol=tol, atol=tol, norm=v -> maximum(abs, v))
near(c::CNDA, re::NDA, im::NDA, tol=1e-13) = near(real(c), re, tol) && near(imag(c), im, tol)
near(c::CNDA, d::CNDA, tol=1e-13) = near(c, real(d), imag(d), tol)

@testset "CNDA" begin
    init!(4, 3, 400)
    # The x1, x2 of test_numeric.cc's CD test, plus a second pair.
    t = davar(1) + 2.0 * davar(2) + 3.0 * davar(3)
    s = 0.5 * davar(1) + 4.0 * davar(2) + 2.7 * davar(3)
    x1, x2, x3, x4 = cos(t), sin(t), cos(s), sin(s)
    a, b = CNDA(x1, x2), CNDA(x3, x4)
    ar, ai, br, bi = x1, x2, x3, x4

    @testset "construction" begin
        @test near(a, x1, x2, 0)
        @test near(CNDA(x1), x1, NDA(0.0), 0)
        @test near(CNDA(1.5 - 2im), NDA(1.5), NDA(-2.0), 0)
        @test near(CNDA(3.0), NDA(3.0), NDA(0.0), 0)
        r = real(a)
        add!(r, r, 1.0)
        @test near(real(a), x1, 0)                        # real returned a copy
        for d in (copy(a), deepcopy(a))
            @test d.ptr != a.ptr && near(d, a, 0)
        end
        @test MiraDAC.env(a) == current_env()
        @test !hasmethod(conj, Tuple{CNDA})                # no conj in C++
    end

    @testset "show" begin
        @test sprint(show, a) == "CNDA(order=4, nvars=3, nonzero=($(nterms(x1)), $(nterms(x2))))"
        full = sprint(show, MIME"text/plain"(), a)
        @test full == MiraDAC.to_string(a)
        lines = split(full, '\n'; keepempty=false)
        @test occursin("V [", lines[1]) && occursin("Base", lines[1])
        @test length(lines) - 2 == count(dense(x1) .!= 0 .|| dense(x2) .!= 0)
    end

    @testset "CNDA operators" begin
        @test near(a + b, ar + br, ai + bi)
        @test near(a - b, ar - br, ai - bi)
        @test near(a * b, ar * br - ai * bi, ar * bi + ai * br)
        n = br * br + bi * bi
        @test near(a / b, (ar * br + ai * bi) / n, (ai * br - ar * bi) / n)
        @test (@inferred a * b) isa CNDA
        @test (@inferred exp(a)) isa CNDA
    end

    @testset "scalar operators" begin
        z = 0.5 - 2.0im
        @test near(a + z, ar + real(z), ai + imag(z))
        @test near(z + a, ar + real(z), ai + imag(z))
        @test near(a - z, ar - real(z), ai - imag(z))
        @test near(z - a, real(z) - ar, imag(z) - ai)
        @test near(a * z, ar * real(z) - ai * imag(z), ar * imag(z) + ai * real(z))
        @test near(z * a, ar * real(z) - ai * imag(z), ar * imag(z) + ai * real(z))
        w = 1 / z
        @test near(a / z, ar * real(w) - ai * imag(w), ar * imag(w) + ai * real(w))
        n = ar * ar + ai * ai
        q = z * CNDA(ar / n, -1.0 * ai / n)
        @test near(z / a, real(q), imag(q))

        @test near(a + 2.0, ar + 2.0, ai)
        @test near(2.0 + a, ar + 2.0, ai)
        @test near(a - 2.0, ar - 2.0, ai)
        @test near(2.0 - a, 2.0 - ar, -1.0 * ai)
        @test near(a * 2.0, ar * 2.0, ai * 2.0)
        @test near(2 * a, ar * 2.0, ai * 2.0)
        @test near(a / 2.0, ar / 2.0, ai / 2.0)
        @test near(2.0 / a, 2.0 * ar / n, -2.0 * ai / n)
        @test near(-a, -ar, -ai, 0)
        @test_throws ArgumentError a / 0.0
    end

    @testset "NDA and complex" begin
        x, z = x1, 0.5 - 2.0im
        @test near(x + z, x + real(z), NDA(imag(z)))
        @test near(z + x, x + real(z), NDA(imag(z)))
        @test near(x - z, x - real(z), NDA(imag(z)))    # as C++ da.h defines it
        @test near(z - x, real(z) - x, NDA(imag(z)))
        @test near(x * z, x * real(z), x * imag(z))
        @test near(z * x, x * real(z), x * imag(z))
        w = 1 / z
        @test near(x / z, x * real(w), x * imag(w))
        @test near(z / x, real(z) / x, imag(z) / x)
        @test x + 2 isa NDA && x * 2.0 isa NDA
    end

    @testset "CNDA and NDA" begin
        x, xc = x1, CNDA(x1)                                # x + 0i
        for op in (+, -, *, /)
            @test norm(real(op(a, x) - op(a, xc))) < 1e-12 && norm(imag(op(a, x) - op(a, xc))) < 1e-12
            @test norm(real(op(x, a) - op(xc, a))) < 1e-12 && norm(imag(op(x, a) - op(xc, a))) < 1e-12
        end
    end

    @testset "in place" begin
        z = 0.5 - 2.0im
        out = CNDA(0.0)
        # The in-place forms keep out's slots: the count does not change. They free queued
        # vectors first, so the garbage is collected and freed before, and none made during, the
        # call.
        function unchanged(f)
            GC.gc()
            MiraDAC.drain!()
            GC.enable(false)
            try
                n = current_env().count
                f()
                return current_env().count == n
            finally
                GC.enable(true)
            end
        end
        for (op, op!) in ((+, add!), (-, sub!), (*, mul!), (/, div!))
            for (l, r) in ((a, b), (a, x3), (x3, a), (a, 2.0), (2.0, a), (a, z), (z, a), (x3, z), (z, x3))
                @test unchanged(() -> @test op!(out, l, r) === out)
                @test near(out, op(l, r), 0)
            end
            c = copy(a)
            op!(c, c, b)                                    # out is an operand
            @test near(c, op(a, b), 0)
        end
        c = 0.3 * a + (0.1 + 0.2im)                         # a's constant part 1 is a branch point
        for f in (sqrt, exp, log, asin, acos, atan, asinh, acosh, atanh)
            @test unchanged(() -> getfield(MiraDAC, Symbol(nameof(f), :!))(out, c))
            @test near(out, f(c), 0)
        end
        for f! in (() -> add!(out, a, b), () -> mul!(out, a, b), () -> mul!(out, a, z),
                   () -> exp!(out, a))
            f!()
            @test (@allocated f!()) == 0
        end
    end

    @testset "pow" begin
        @test near(a^3, a * a * a)
        r = a^0.5
        @test near(r * r, a, 1e-12)
    end

    @testset "functions" begin
        e = exp(a)
        @test near(e, exp(ar) * cos(ai), exp(ar) * sin(ai))
        r = sqrt(a)
        @test near(r * r, a)
        @test near(exp(log(a)), a)
        @test near(a^0.5, r)
        @test abs(a) == max(norm(ar), norm(ai))
        c = 0.3 * a + (0.1 + 0.2im)                         # away from asin's branch point at 1
        @test near(asin(c) + acos(c), NDA(pi / 2), NDA(0.0), 1e-12)
        @test_throws MethodError sin(a)                     # no complex sin in C++
        for (f, c0) in ((atan, 0.3), (asinh, 0.3), (atanh, 0.3), (acosh, 1.7), (asin, 0.3), (acos, 0.3))
            x = c0 + davar(1) + 0.5 * davar(2)
            v = f(CNDA(x))
            @test near(v, f(x), NDA(0.0), 1e-12)
            @test con(real(v)) ≈ real(f(complex(c0)))
        end
    end

    @testset "env mismatch" begin
        clear!()
        init!(4, 3, 400)
        @test_throws EnvError a + 1.0
        @test_throws EnvError CNDA(davar(1)) + CNDA(1.0) + a
        @test_throws EnvError CNDA(davar(1), x1)
        @test sprint(show, a) == "CNDA(cleared env)"
    end

    init!(4, 3, 400)
    t = davar(1) + 2.0 * davar(2) + 3.0 * davar(3)
    s = 0.5 * davar(1) + 4.0 * davar(2) + 2.7 * davar(3)
    x1, x2, x3, x4 = cos(t), sin(t), cos(s), sin(s)
    y1, y2 = CNDA(x1, x2), CNDA(x3, x4)

    @testset "cd_calculation files" begin
        for (i, op) in enumerate((+, -, *, /))
            @test compare_cd_with_file("cd_calculation_$(i-1).txt", op(y1, y2), 1e-14)
        end
        @test !compare_cd_with_file("cd_calculation_2.txt", y1 + y2, 1e-14)
    end

    @testset "CNDAList" begin
        l = CNDAList([y1])
        push!(l, y2)
        @test length(l) == 2 && l isa AbstractVector{CNDA}
        l[1] = l[1] + 1.0
        @test con(real(l[1])) ≈ con(x1) + 1.0
        n = current_env().count
        l[2] = y1
        @test current_env().count == n                      # set copies into the element's slots
        @test near(l[2], y1, 0)
        @test_throws BoundsError l[3]
        @test length(CNDAList()) == 0
        @test_throws EnvError MiraDAC.env(CNDAList())
    end

    @testset "cd_composition files" begin
        for as in (identity, collect)
            mmap = NDAList([x1, x2])
            cnmap = CNDAList([y1, y2, y1 * y2])
            cmmap = CNDAList([CNDA(x1, exp(x1)), CNDA(x2, exp(x2))])
            out = compose(as(mmap), as(cnmap))
            @test out isa CNDAList && length(out) == 2
            @test all(compare_cd_with_file("da_composition_cd_$(i-1).txt", out[i], 1e-14) for i in 1:2)
            out = compose(as(cmmap), as(cnmap))
            @test all(compare_cd_with_file("cd_composition_cd_$(i-1).txt", out[i], 1e-14) for i in 1:2)
            push!(mmap, x1 + 0.33 * x2)
            out = compose(as(cmmap), as(mmap))
            @test all(compare_cd_with_file("cd_composition_da_$(i-1).txt", out[i], 1e-14) for i in 1:2)
        end
        @test_throws ArgumentError compose([x1, x2], [y1, y2])
        @test_throws ArgumentError compose([y1], [x1])
        @test length(compose(CNDA[], [y1, y2, y1])) == 0
    end

    @testset "compose at complex points" begin
        mmap = [x1, x2]
        pt = [0.1 + 0.2im, -0.3 + 0.05im, 0.2 - 0.1im]
        vals = compose(mmap, pt)
        out = compose(mmap, [CNDA(z) for z in pt])
        for (v, c) in zip(vals, out)
            @test v ≈ complex(con(real(c)), con(imag(c))) rtol = 1e-14
        end
    end
end
