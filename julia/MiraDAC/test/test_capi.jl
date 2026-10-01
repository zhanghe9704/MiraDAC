# T2.1: status checking and the allocation retry.
@testset "capi" begin
    init!(4, 3, 400)

    @testset "check" begin
        @test MiraDAC.check(MiraDAC.MDAC_OK) === nothing
        @test_throws EnvError MiraDAC.check(MiraDAC.MDAC_ERR_ENV)
        @test_throws ArgumentError MiraDAC.check(MiraDAC.MDAC_ERR_VALUE)
        @test_throws DomainError MiraDAC.check(MiraDAC.MDAC_ERR_VALUE, 2.0)
        @test_throws BoundsError MiraDAC.check(MiraDAC.MDAC_ERR_INDEX)
        @test_throws PoolExhaustedError MiraDAC.check(MiraDAC.MDAC_ERR_POOL)
        @test_throws ErrorException MiraDAC.check(MiraDAC.MDAC_ERR_RUNTIME)
        @test_throws ErrorException MiraDAC.check(MiraDAC.MDAC_ERR_UNSUPPORTED)
        @test EnvError <: MiraDACError && PoolExhaustedError <: MiraDACError
        # The message is mdac_last_error().
        err = try set_eps!(-1.0) catch e; e end
        @test err isa ArgumentError && err.msg == "eps must be positive"
    end

    @testset "alloc_call retry" begin
        n0 = MiraDAC.POOL_RETRIES[]
        calls = Ref(0)
        f = () -> (calls[] += 1; calls[] == 1 ? MiraDAC.MDAC_ERR_POOL : MiraDAC.MDAC_OK)
        @test MiraDAC.alloc_call(f) === nothing
        @test calls[] == 2
        @test MiraDAC.POOL_RETRIES[] == n0 + 1

        calls[] = 0
        @test_throws PoolExhaustedError MiraDAC.alloc_call(() -> (calls[] += 1; MiraDAC.MDAC_ERR_POOL))
        @test calls[] == 1 + MiraDAC.retry_rounds()
        @test_throws ArgumentError MiraDAC.alloc_call(() -> MiraDAC.MDAC_ERR_VALUE)
    end
end
