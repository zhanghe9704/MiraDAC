# T2.5: the pool, the GC and finalizers (A.5).

# DAEnv(order, nvars, poolsize) is Stage 7; until then, the C API directly.
function make_env(order, nvars, poolsize)
    r = Ref{MiraDAC.Handle}(C_NULL)
    MiraDAC.check(MiraDAC.mdac_env_make(order, nvars, poolsize, 0, r))
    return DAEnv(r[])
end

@testset "memory" begin
    @testset "retry recovers" begin
        init!(4, 3, 64)
        n0 = MiraDAC.POOL_RETRIES[]
        x = NDA(0.0)
        for _ in 1:100_000
            x = x * 0.5 + davar(1)
        end
        @test con(x) == 0.0 && coeff(x, [1, 0, 0]) ≈ 2.0
        @test MiraDAC.POOL_RETRIES[] > n0
    end

    # The C++ engine frees its temporaries when a kernel throws (commit a6be310), so a retry
    # inside exp or compose recovers too.
    @testset "retry recovers in functions with temporaries" begin
        init!(4, 3, 64)
        n0 = MiraDAC.POOL_RETRIES[]
        a = 0.3 + davar(1)
        for _ in 1:10_000
            exp(a)
        end
        @test con(exp(a)) ≈ exp(0.3)
        m = [davar(1) + 0.1, davar(2) * davar(3), davar(3)]
        n = [davar(2) + 0.2, davar(1) - 1.0, davar(3) * 2.0]
        local r
        for _ in 1:10_000
            r = compose(m, n)
        end
        @test con(r[1]) ≈ 0.3 && coeff(r[2], [1, 0, 1]) ≈ 2.0
        @test MiraDAC.POOL_RETRIES[] > n0
    end

    @testset "exhaustion throws" begin
        init!(4, 3, 64)
        live = NDA[]
        @test_throws PoolExhaustedError for _ in 1:100
            push!(live, NDA(1.0))
        end
        @test length(live) < 64
        empty!(live)
    end

    @testset "use after clear" begin
        init!(4, 3, 64)
        v = davar(1)
        w = [davar(2) for _ in 1:10]
        clear!()
        @test_throws EnvError v + v
        @test_throws EnvError con(v)
        @test_throws EnvError add!(v, v, v)
        v = w = nothing
        GC.gc(); GC.gc()
        MiraDAC.drain!()
        init!(4, 3, 400)
        @test con(davar(1) + 1.0) == 1.0
    end

    @testset "finalizers on other threads ($(Threads.nthreads()) threads)" begin
        envs = [make_env(4, 3, 64) for _ in 1:4]
        base = [e.count for e in envs]
        tasks = map(envs) do e
            Threads.@spawn begin
                prev = Ref{MiraDAC.Handle}(C_NULL)
                MiraDAC.check(MiraDAC.mdac_env_exchange(e, prev))
                try
                    x = NDA(0.0)
                    for _ in 1:20_000
                        x = x * 0.5 + davar(1)
                    end
                    coeff(x, [1, 0, 0])
                finally
                    MiraDAC.mdac_env_exchange(prev[], Ref{MiraDAC.Handle}(C_NULL))
                end
            end
        end
        for _ in 1:3  # collections from the main thread while the tasks run
            GC.gc()
            sleep(0.01)
        end
        @test all(t -> fetch(t) ≈ 2.0, tasks)
        GC.gc(); GC.gc()
        MiraDAC.drain!()
        @test [e.count for e in envs] == base
        foreach(e -> MiraDAC.check(MiraDAC.mdac_env_close(e)), envs)
    end
end
