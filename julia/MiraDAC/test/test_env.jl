# T2.4: the default env, the order and eps.
@testset "env" begin
    e = init!(4, 3, 400)
    @test e == default_env() == current_env()
    @test (e.order, e.max_order, e.nvars, e.full_length, e.poolsize) == (4, 4, 3, 35, 400)
    @test !e.retired

    @testset "with_order" begin
        r = with_order(3) do
            @test current_env().order == 3
            with_order(2) do
                @test current_env().order == 2
                @test nterms(davar(1)^3) == 0
            end
            @test current_env().order == 3
            :done
        end
        @test r == :done
        @test current_env().order == 4
        with_order(3) do
            @test_throws ErrorException with_order(1) do
                error("boom")
            end
            @test current_env().order == 3
        end
        @test current_env().order == 4
        @test_throws ArgumentError with_order(() -> nothing, 5)
        @test current_env().order == 4
    end

    @testset "eps" begin
        old = get_eps()
        try
            set_eps!(1e-10)
            @test get_eps() == 1e-10
            @test_throws ArgumentError set_eps!(0.0)
        finally
            set_eps!(old)
        end
    end

    @testset "clear" begin
        v = davar(1)
        clear!()
        @test e.retired
        @test_throws EnvError current_env()
        @test_throws EnvError default_env()
        @test_throws EnvError NDA(1.0)
        @test sprint(show, v) == "NDA(cleared env)"
        init!(4, 3, 400)
    end
end
