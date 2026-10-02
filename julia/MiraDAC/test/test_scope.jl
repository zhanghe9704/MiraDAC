# A.5b: dascope frees temporaries at scope exit; keep! and returned values survive.

@testset "dascope" begin
    init!(4, 3, 64)
    e = current_env()

    @testset "returned values survive, temporaries are freed" begin
        x = davar(1) + 0.5
        n0 = e.count
        y = dascope() do
            t = exp(x) * sin(x)          # temporaries
            t + 1.0
        end
        @test e.count == n0 + 1          # only y is left; no GC needed
        @test con(y) ≈ exp(0.5) * sin(0.5) + 1.0
        @test coeff(y, [1, 0, 0]) ≈ exp(0.5) * (sin(0.5) + cos(0.5))
        @test dascope(() -> 42) == 42
        @test dascope(() -> con(exp(x))) ≈ exp(0.5)
        @test e.count == n0 + 1
    end

    @testset "nested structures returned survive" begin
        n0 = e.count
        r = dascope() do
            a, b = davar(1), davar(2)
            junk = a * b * 3.0
            (pair=(a + 1.0, [b * 2.0, a * b]), d=Dict(:c => CNDA(a, b)), l=NDAList([a, b]),
             any=Any[Any[a - b]], n=3)
        end
        @test e.count == n0 + 8     # a+1, b*2, a*b, CNDA (2), list (2), a-b
        @test con(r.pair[1]) == 1.0
        @test coeff(r.pair[2][1], [0, 1, 0]) == 2.0
        @test coeff(r.pair[2][2], [1, 1, 0]) == 1.0
        @test coeff(real(r.d[:c]), [1, 0, 0]) == 1.0
        @test coeff(r.l[2], [0, 1, 0]) == 1.0
        @test coeff(r.any[1][1], [0, 1, 0]) == -1.0
        @test r.n == 3
        # A partly filled array (#undef slots) is a valid value.
        nu = e.count
        out = dascope() do
            o = Vector{NDA}(undef, 3)
            o[1] = davar(1) + 1.0
            o
        end
        @test e.count == nu + 1 && !isassigned(out, 2)
        @test con(out[1]) == 1.0
        # Lazy wrappers keep what they wrap, without indexing it (no copies, no adjoint(::NDA)).
        nv = e.count
        vl = dascope(() -> (l = NDAList([davar(1), davar(2)]); view(l, 1:2)))
        @test e.count == nv + 2          # the list only: no element copies during the keep
        @test coeff(vl[2], [0, 1, 0]) == 1.0
        va = dascope(() -> [davar(1), davar(2) * 2.0]')
        @test va isa LinearAlgebra.Adjoint && coeff(parent(va)[2], [0, 1, 0]) == 2.0
        vr = dascope(() -> reshape(view([davar(1), davar(2), davar(3), davar(1)], :), 2, 2))
        @test coeff(vr[2, 2], [1, 0, 0]) == 1.0
        # At the outermost scope, kept objects are left to the GC (a function: no global roots).
        dropped() = (dascope(() -> (davar(1) + 1.0, [davar(2) * 2.0])); nothing)
        n1 = e.count
        dropped()
        GC.gc(); GC.gc(); MiraDAC.drain!()
        @test e.count == n1
    end

    @testset "keep!" begin
        n0 = e.count
        local k
        v = dascope() do
            k = keep!(davar(2) * 3.0)
            w = davar(3)
            keep!((w,))
            nothing
        end
        @test v === nothing
        @test e.count == n0 + 2
        @test coeff(k, [0, 1, 0]) == 3.0
        @test keep!(5) == 5                    # outside a scope: nothing to do
    end

    @testset "freed objects throw FreedObjectError" begin
        local t, l, c
        dascope() do
            t = davar(1) * 2.0
            l = NDAList([t])
            c = CNDA(t)
            nothing
        end
        @test t.ptr == C_NULL
        @test_throws FreedObjectError t + 1.0
        @test_throws FreedObjectError con(t)
        @test_throws FreedObjectError add!(t, davar(1), davar(2))
        @test_throws FreedObjectError NDAList([t])
        @test_throws FreedObjectError length(l)
        @test_throws FreedObjectError c * 2.0
        @test occursin("keep!", sprint(showerror, try con(t) catch err; err end))
    end

    @testset "nesting moves kept objects outward" begin
        n0 = e.count
        local inner, tmp
        outer = dascope() do
            inner = dascope() do
                tmp = davar(1) * 5.0
                tmp + davar(2)
            end
            @test e.count == n0 + 1               # inner temporaries freed
            @test tmp.ptr == C_NULL
            @test coeff(inner, [0, 1, 0]) == 1.0  # inner result still alive
            inner * 2.0
        end
        @test inner.ptr == C_NULL                 # freed by the outer scope
        @test coeff(outer, [1, 0, 0]) == 10.0
        @test e.count == n0 + 1
    end

    @testset "an exception frees temporaries and rethrows" begin
        n0 = e.count
        local t, k
        @test_throws ErrorException dascope() do
            t = davar(1) + 1.0
            k = keep!(davar(2) + 1.0)
            error("boom")
        end
        @test e.count == n0 + 1
        @test t.ptr == C_NULL
        @test con(k) == 1.0
        @test_throws DomainError dascope(() -> asin(NDA(2.0)))
        @test e.count == n0 + 1
        # The scope stack is unwound: new objects get finalizers again.
        @test MiraDAC.current_scope() === nothing
        @test MiraDAC.ACTIVE_SCOPES[] == 0
    end

    @testset "scopes are per task" begin
        dascope() do
            t = fetch(Threads.@spawn MiraDAC.current_scope())
            @test t === nothing
            @test MiraDAC.current_scope() !== nothing
        end
    end

    @testset "100 000 operations in scopes, pool 64: no retry" begin
        GC.gc(); GC.gc(); MiraDAC.drain!()
        init!(4, 3, 64)
        e = current_env()
        base = e.count
        n0 = MiraDAC.POOL_RETRIES[]
        x = NDA(0.0)
        s = 0.0
        for _ in 1:50_000
            dascope() do
                add!(x, x * 0.5, davar(1))           # temporaries, result in place
            end
            s += dascope(() -> con(exp(x * 0.5 + 0.1)))
        end
        @test coeff(x, [1, 0, 0]) ≈ 2.0
        @test s ≈ 50_000 * exp(0.1)
        @test MiraDAC.POOL_RETRIES[] == n0
        @test e.count == base + 1                   # x
        # One scope per batch of 16 operations (3 objects each: 48 live at most).
        y = dascope() do
            for _ in 1:6_250
                dascope() do
                    for _ in 1:16
                        x * 0.5 + davar(1)
                    end
                end
            end
            x * 1.0
        end
        @test MiraDAC.POOL_RETRIES[] == n0
        @test e.count == base + 2
        @test coeff(y, [1, 0, 0]) ≈ 2.0
    end

    @testset "errors raised in a scope keep their message" begin
        init!(4, 3, 64)
        err = try
            dascope(() -> asin(NDA(2.0)))
            nothing
        catch e
            e
        end
        @test err isa DomainError
        msg = sprint(showerror, err)          # used to throw FreedObjectError
        @test occursin("asin", msg)
        @test occursin("freed", msg)
        v = Ref{Any}()
        dascope() do
            v[] = NDA(1.0)
            nothing
        end
        @test sprint(show, v[]) == "NDA(freed)"
        @test sprint(show, MIME"text/plain"(), v[]) == "NDA(freed)"
    end

    if MiraDAC.HAS_SYMBOLIC
        @testset "SDA and CSDA are scoped too" begin
            init!(3, 2, 64)
            a, = dasymbols("a")
            tmp = Ref{Any}()
            ctmp = Ref{Any}()
            y = dascope() do
                t = a * sdavar(1) + SDA(1.0)
                tmp[] = t
                ctmp[] = CSDA(t, t)
                exp(t)
            end
            @test sprint(show, tmp[]) == "SDA(freed)"
            @test sprint(show, ctmp[]) == "CSDA(freed)"
            @test_throws FreedObjectError tmp[] + 1.0
            @test con(evaluate(y, Dict(a => 0.0))) ≈ exp(1.0)
        end
    end
end
