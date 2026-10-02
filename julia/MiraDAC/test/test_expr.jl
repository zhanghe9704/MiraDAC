# T5.3: SymExpr, ported from python/tests/test_expr.py (without the sympy tests).

@testset "SymExpr" begin
    @testset "constructors and show" begin
        @test string(SymExpr(3)) == "3"
        @test Float64(SymExpr(0.25)) == 0.25
        @test sprint(show, SymExpr("a + 2*b")) == "SymExpr(\"a + 2*b\")"
        @test SymExpr(big(2)^100) == SymExpr(string(big(2)^100))
        @test SymExpr(typemin(Int64)) == SymExpr(string(typemin(Int64)))
        @test SymExpr(3) == 3 && SymExpr(1.5) == 1.5
        @test_throws ArgumentError SymExpr("a +* b")
        a = SymExpr("a")
        for c in (copy(a), deepcopy(a))
            @test c.ptr != a.ptr && c == a
        end
    end

    @testset "dasymbols" begin
        a, b, c = dasymbols("a b c")
        @test (string(a), string(b), string(c)) == ("a", "b", "c")
        @test dasymbols("x, y") == dasymbols("x y")
        @test length(dasymbols("  z ")) == 1
    end

    @testset "arithmetic" begin
        a, b = dasymbols("a b")
        @test a + b == SymExpr("a + b")
        @test 2 * a - b / 3 == SymExpr("2*a - b/3")
        @test 1 + a == a + 1
        @test 1.5 * a == SymExpr("1.5*a")
        @test 1 / a == a^-1
        @test a^b == SymExpr("a**b")
        @test 2^a == SymExpr("2**a")
        @test a^0.5 == SymExpr("a**0.5")
        @test -a == SymExpr("-a")
        @test +a == a
        @test a - a == 0
        @test a - 0.5 == SymExpr("-0.5 + a") && 0.5 - a == SymExpr("0.5 - a")
    end

    @testset "Float64, ==, hash" begin
        a, b = dasymbols("a b")
        @test Float64(SymExpr("sin(1) + 2")) ≈ sin(1) + 2
        @test_throws ArgumentError Float64(a + 1)
        @test a + b == b + a
        @test a != b
        @test hash(a + b) == hash(b + a)
        @test length(Set([a, b, a + 0])) == 2
        @test (a == "a") === false
    end

    @testset "methods" begin
        a, b = dasymbols("a b")
        e = (a + b)^2
        @test MiraDAC.expand(e) == a^2 + 2 * a * b + b^2
        @test MiraDAC.subs(e, Dict(a => 1.0, b => 2)) == SymExpr(9.0)
        @test MiraDAC.subs(e, Dict(a => b)) == 4 * b^2
        @test diff(e, a) == 2 * (a + b)
        @test_throws ArgumentError diff(e, a + b)
        @test MiraDAC.free_symbols(e) == Set([a, b])
        @test MiraDAC.free_symbols(SymExpr(3)) == Set{SymExpr}()
        @test iszero(e - MiraDAC.expand(e))
        @test !iszero(e)
        @test simplify(SymExpr("sin(a)**2 + cos(a)**2")) isa SymExpr
        @test simplify(a * 2 / 2) == a
    end
end
