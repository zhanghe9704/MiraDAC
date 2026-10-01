# T2.2, T2.3: NDA construction, inspection, operators and functions.

# Port of da::compare_da_with_file (src/base.cpp): read the reference table, then every
# coefficient of (ref - v) ./ v (by v's coefficient where it is non-zero) is below eps.
function read_da_file(path)
    terms = Dict{Vector{Int},Float64}()
    reading = false
    for line in eachline(path)
        line = strip(line)
        isempty(line) && continue
        if reading
            w = split(line)
            terms[parse.(Int, w[3:end-1])] = parse(Float64, w[2])
        elseif isempty(replace(line, '-' => ""))
            reading = true
        end
    end
    return terms
end

function nda_terms(v::NDA)
    terms = Dict{Vector{Int},Float64}()
    exps = Vector{Cint}(undef, MiraDAC.env(v).nvars)
    x = Ref{Cdouble}()
    for i in 0:length(coeffs(v))-1
        MiraDAC.check(MiraDAC.mdac_nda_index_term(v, i, exps, x))
        x[] == 0 || (terms[Int.(exps)] = x[])
    end
    return terms
end

ref_path(name) = joinpath(@__DIR__, "..", "..", "..", "test", name)

# da::compare_da_vectors on term tables.
function compare_terms(ref, got, eps)
    return all(union(keys(ref), keys(got))) do k
        r = get(ref, k, 0.0) - get(got, k, 0.0)
        b = get(got, k, 0.0)
        abs(abs(b) > floatmin(Float64) ? r / b : r) < eps
    end
end

compare_da_with_file(name, v::NDA, eps) = compare_terms(read_da_file(ref_path(name)), nda_terms(v), eps)

@testset "NDA" begin
    init!(4, 3, 400)
    x, y, z = davar(1), davar(2), davar(3)

    @testset "construction and inspection" begin
        v = 1.0 + x + 2.0 * y + 5.0 * z
        @test coeffs(NDA(coeffs(v))) == coeffs(v)
        @test coeffs(NDA([1, 2, 3])) == [1.0, 2.0, 3.0]
        @test con(NDA(2.5)) == 2.5
        @test coeff(v, [0, 1, 0]) == 2.0
        @test coeff(v, (0, 0, 1)) == 5.0
        @test coeff(v, [0]) == 1.0
        @test nterms(v) == 4
        @test norm(v) == 5.0
        @test iszero(NDA(0.0))
        @test !iszero(x)
        @test iszero(NDA(1e-5); eps=1e-4)
        @test MiraDAC.env(v) == current_env()
        @test_throws BoundsError davar(0)
        @test_throws BoundsError davar(4)
        @test_throws ArgumentError coeff(v, [-1, 0, 0])

        w = copy(v)
        @test w.ptr != v.ptr && coeffs(w) == coeffs(v)
        d = deepcopy([v, v])
        @test d[1] === d[2] && d[1].ptr != v.ptr && coeffs(d[1]) == coeffs(v)

        @test sprint(show, x) == "NDA(order=4, nvars=3, nonzero=1)"
        full = sprint(show, MIME"text/plain"(), x)
        @test full == MiraDAC.to_string(x)
        @test occursin(r"1\s+1\.000000000000000e\+00\s+1 0 0\s+1", full)
    end

    a, b = 1.0 + x, 2.0 + y
    c(v, e...) = coeff(v, collect(e))

    @testset "operators" begin
        s = a + b
        @test (con(s), c(s, 1, 0, 0), c(s, 0, 1, 0)) == (3.0, 1.0, 1.0)
        s = a - b
        @test (con(s), c(s, 1, 0, 0), c(s, 0, 1, 0)) == (-1.0, 1.0, -1.0)
        s = a * b
        @test (con(s), c(s, 1, 0, 0), c(s, 0, 1, 0), c(s, 1, 1, 0), nterms(s)) == (2.0, 2.0, 1.0, 1.0, 4)
        s = a / b
        @test (con(s), c(s, 1, 0, 0), c(s, 0, 1, 0), c(s, 1, 1, 0), c(s, 0, 2, 0)) ==
              (0.5, 0.5, -0.25, -0.25, 0.125)
        @test coeffs(a + 2.5) == coeffs(2.5 + a) == [3.5, 1.0]
        @test coeffs(a - 2.5) == [-1.5, 1.0]
        @test coeffs(2.5 - a) == [1.5, -1.0]
        @test coeffs(a * 3) == coeffs(3 * a) == [3.0, 3.0]
        @test coeffs(a / 2) == [0.5, 0.5]
        @test [c(2 / a, k, 0, 0) for k in 0:4] == [2.0, -2.0, 2.0, -2.0, 2.0]
        @test coeffs(-a) == [-1.0, -1.0]
        @test [c(a^2, k, 0, 0) for k in 0:2] == [1.0, 2.0, 1.0]
        @test [c(a^0.5, k, 0, 0) for k in 0:4] ≈ [1, 1/2, -1/8, 1/16, -5/128]
        @test_throws ArgumentError a / 0.0
        @test (@inferred a + b) isa NDA
        @test (@inferred a * 2.0) isa NDA
        @test (@inferred exp(a)) isa NDA
    end

    @testset "functions" begin
        # f(c + x): constant f(c), first-order coefficient f'(c).
        cases = [(sqrt, 0.3, 1 / (2sqrt(0.3))), (exp, 0.3, exp(0.3)), (log, 0.3, 1 / 0.3),
                 (sin, 0.3, cos(0.3)), (cos, 0.3, -sin(0.3)), (tan, 0.3, 1 + tan(0.3)^2),
                 (asin, 0.3, 1 / sqrt(1 - 0.09)), (acos, 0.3, -1 / sqrt(1 - 0.09)),
                 (atan, 0.3, 1 / 1.09), (sinh, 0.3, cosh(0.3)), (cosh, 0.3, sinh(0.3)),
                 (tanh, 0.3, 1 - tanh(0.3)^2), (asinh, 0.3, 1 / sqrt(1.09)),
                 (acosh, 1.5, 1 / sqrt(1.25)), (atanh, 0.3, 1 / (1 - 0.09)),
                 (erf, 0.3, 2 / sqrt(pi) * exp(-0.09))]
        erf_03 = 0.328626759459127427638914300867  # erf(0.3)
        for (f, x0, d) in cases
            v = f(x0 + x)
            @test con(v) ≈ (f === erf ? erf_03 : f(x0))
            @test c(v, 1, 0, 0) ≈ d
            out = NDA(0.0)
            @test getfield(MiraDAC, Symbol(nameof(f), :!))(out, x0 + x) === out
            @test coeffs(out) == coeffs(v)
        end
        @test abs(1.0 - 3.0 * y) == 3.0
        @test_throws DomainError asin(NDA(2.0))
        @test_throws DomainError asin!(NDA(0.0), NDA(2.0))
    end

    @testset "in place" begin
        out = NDA(0.0)
        @test coeffs(add!(out, a, b)) == coeffs(a + b)
        @test coeffs(sub!(out, a, b)) == coeffs(a - b)
        @test coeffs(mul!(out, a, b)) == coeffs(a * b)
        @test coeffs(div!(out, a, b)) == coeffs(a / b)
        @test coeffs(add!(out, a, 2.0)) == coeffs(a + 2.0)
        @test coeffs(sub!(out, 2.0, a)) == coeffs(2.0 - a)
        @test coeffs(mul!(out, a, 3.0)) == coeffs(a * 3.0)
        @test coeffs(div!(out, 2.0, a)) == coeffs(2.0 / a)
        n = current_env().count
        add!(out, a, b); mul!(out, a, b); exp!(out, a)
        @test current_env().count == n
        for f! in (() -> add!(out, a, b), () -> mul!(out, a, b), () -> add!(out, a, 2.0),
                   () -> exp!(out, a))
            f!()
            @test (@allocated f!()) == 0
        end
        @test (add!(out, a, b); @allocated add!(out, a, b)) == 0
    end

    @testset "reference files" begin
        v = 1.0 + x + 2.0 * y + 5.0 * z
        @test compare_da_with_file("sqrt_da.txt", sqrt(v), 1e-14)
        @test compare_da_with_file("log_da.txt", log(v), 1e-14)
        @test compare_da_with_file("exp_da.txt", exp(v), 1e-14)
        @test compare_da_with_file("pow0p3_da.txt", v^0.3, 1e-14)
        @test compare_da_with_file("pow3_da.txt", v^3, 1e-14)
        @test !compare_da_with_file("exp_da.txt", exp(v) * 1.001, 1e-14)
    end
end
